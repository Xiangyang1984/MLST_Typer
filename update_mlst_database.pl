#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename);
use JSON::PP qw(decode_json encode_json);
use MLSTTools qw(load_scheme_table select_scheme list_species list_schemes safe_name which_cmd current_time);

my %opt = (
    scheme_table => "$FindBin::Bin/data/mlst_schemes_all.tab",
    db_root      => "$FindBin::Bin/database",
    threads      => 4,
    auth_dir     => "$FindBin::Bin/auth",
);
GetOptions(
    'species=s'          => \$opt{species},
    'scheme=s'           => \$opt{scheme},
    'scheme-description=s'=>\$opt{description},
    'scheme-table=s'     => \$opt{scheme_table},
    'db-root=s'          => \$opt{db_root},
    'threads=i'          => \$opt{threads},
    'auth-dir=s'         => \$opt{auth_dir},
    'list-species'       => \$opt{list_species},
    'list-schemes'       => \$opt{list_schemes},
    'force'              => \$opt{force},
    'help|h'             => \$opt{help},
) or die usage();

if ($opt{help}) { print usage(); exit 0; }
my $rows = load_scheme_table($opt{scheme_table});
if ($opt{list_species}) {
    print "Species\tNumber_of_schemes\n";
    print $_->[0], "\t", $_->[1], "\n" for @{list_species($rows)};
    exit 0;
}
if ($opt{list_schemes}) {
    die "--species is required with --list-schemes\n" unless $opt{species};
    print "database\tspecies\tscheme_description\tscheme\tURI\n";
    my $m = list_schemes($rows,$opt{species});
    print join("\t", @$_{qw(database species scheme_description scheme URI)}), "\n" for @$m;
    exit 0;
}

my $sel = select_scheme(rows=>$rows, species=>$opt{species}, scheme=>$opt{scheme}, description=>$opt{description}, interactive=>1);
my $scheme = $sel->{scheme};
my $dest = "$opt{db_root}/$scheme";
my $stage = "$opt{db_root}/.${scheme}.update.$$";
make_path($opt{db_root});
remove_tree($stage) if -d $stage;
make_path("$stage/loci");

print "[", current_time(), "] Selected MLST scheme\n";
print "  Species: $sel->{species}\n  Scheme: $scheme\n  Description: $sel->{scheme_description}\n  Source: $sel->{database}\n  URI: $sel->{URI}\n";
print "  OAuth tokens: $opt{auth_dir} (one global token set for all PubMLST schemes)\n" if $sel->{URI} =~ m{^https?://rest\.pubmlst\.org/}i;
print "[", current_time(), "] Reading scheme metadata...\n";

my $scheme_json_file = "$stage/scheme.json";
fetch_url($sel->{URI}, $scheme_json_file);
my $scheme_obj = decode_json(slurp($scheme_json_file));
my @locus_uris = ref($scheme_obj->{loci}) eq 'ARRAY' ? @{$scheme_obj->{loci}} : ();
if (!@locus_uris) {
    my $lf = "$stage/loci.json";
    fetch_url($sel->{URI} . '/loci', $lf);
    my $lo = decode_json(slurp($lf));
    @locus_uris = @{$lo->{loci} || []};
}
die "No loci were returned for scheme $scheme.\n" unless @locus_uris;

my @loci = map { locus_from_uri($_) } @locus_uris;
open my $lout, '>', "$stage/loci.txt" or die $!;
print $lout "$_\n" for @loci;
close $lout;

my $profiles_url = $scheme_obj->{profiles_csv} || ($sel->{URI} . '/profiles_csv');
print "[", current_time(), "] Downloading profiles...\n";
fetch_url($profiles_url, "$stage/profiles.tsv");
validate_profiles("$stage/profiles.tsv", \@loci);

print "[", current_time(), "] Downloading ", scalar(@loci), " loci (max $opt{threads} parallel requests)...\n";
my @jobs;
for my $i (0..$#loci) {
    push @jobs, [$loci[$i], $locus_uris[$i]];
}
parallel_download(\@jobs, $stage, $opt{threads});

print "[", current_time(), "] Building canonical allele FASTA and index...\n";
open my $all, '>', "$stage/alleles.fasta" or die $!;
open my $idx, '>', "$stage/allele_index.tsv" or die $!;
print $idx "blast_id\tlocus\tallele\toriginal_id\n";
my $counter=0;
for my $locus (@loci) {
    my $safe = safe_name($locus);
    my $f = "$stage/loci/$safe.fas";
    open my $fh, '<', $f or die "Cannot open $f: $!\n";
    my ($id,$seq)=('','');
    while (<$fh>) {
        chomp;
        if (/^>(\S+)/) {
            emit_allele($all,$idx,\$counter,$locus,$id,$seq) if $id ne '';
            ($id,$seq)=($1,'');
        } else { s/\s+//g; $seq .= $_; }
    }
    emit_allele($all,$idx,\$counter,$locus,$id,$seq) if $id ne '';
    close $fh;
}
close $all; close $idx;
die "No allele sequences were downloaded.\n" unless -s "$stage/alleles.fasta";

open my $info, '>', "$stage/scheme.info.tsv" or die $!;
print $info "key\tvalue\n";
for my $p (
 ['database',$sel->{database}], ['species',$sel->{species}], ['scheme_description',$sel->{scheme_description}],
 ['scheme',$scheme], ['URI',$sel->{URI}], ['locus_count',scalar(@loci)], ['allele_count',$counter], ['updated_at',current_time()]
) { print $info "$p->[0]\t$p->[1]\n"; }
close $info;

# Atomic-ish replacement: only replace the previous scheme after a complete staged update.
my $backup = "$opt{db_root}/.${scheme}.backup.$$";
if (-d $dest) { rename($dest,$backup) or die "Cannot backup existing $dest: $!\n"; }
rename($stage,$dest) or do {
    rename($backup,$dest) if -d $backup;
    die "Cannot install updated database to $dest: $!\n";
};
remove_tree($backup) if -d $backup;
open my $log, '>>', "$opt{db_root}/Update_time.txt" or die $!;
print $log current_time(), "\t$sel->{species}\t$scheme\t$sel->{scheme_description}\n";
close $log;
print "[", current_time(), "] MLST database update completed: $dest\n";

sub emit_allele {
    my ($all,$idx,$counter_ref,$locus,$id,$seq)=@_;
    $seq =~ s/[^A-Za-z]//g; $seq=uc($seq);
    return if $seq eq '';
    $$counter_ref++;
    my $bid = sprintf('A%09d',$$counter_ref);
    my $allele = extract_allele_id($locus,$id);
    print $all ">$bid\n";
    for (my $i=0;$i<length($seq);$i+=60) { print $all substr($seq,$i,60),"\n"; }
    print $idx join("\t",$bid,$locus,$allele,$id),"\n";
}

sub extract_allele_id {
    my ($locus,$id)=@_;
    my $a=$id;
    if ($a =~ /^\Q$locus\E[_-](.+)$/i) { return $1; }
    if ($a =~ /[_-]([^_-]+)$/ && $1 =~ /^[A-Za-z0-9.]+$/) { return $1; }
    return $a;
}

sub locus_from_uri {
    my ($u)=@_;
    $u =~ s/[?#].*$//;
    $u =~ s{/$}{};
    my $x = basename($u);
    $x =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/eg;
    return $x;
}

sub parallel_download {
    my ($jobs,$stage,$threads)=@_;
    $threads = 1 if !$threads || $threads < 1;
    $threads = 4 if $threads > 4; # PubMLST asks clients not to exceed 4 simultaneous requests.
    my @queue=@$jobs;
    my %children;
    my $done=0; my $total=@queue;
    while (@queue || %children) {
        while (@queue && keys(%children) < $threads) {
            my $job=shift @queue;
            my ($locus,$uri)=@$job;
            my $pid=fork();
            die "fork failed: $!\n" unless defined $pid;
            if ($pid==0) {
                my $safe=safe_name($locus);
                eval { fetch_url($uri . '/alleles_fasta', "$stage/loci/$safe.fas"); 1 } or do { warn $@; exit 2; };
                exit 0;
            }
            $children{$pid}=$locus;
        }
        my $pid=wait();
        last if $pid < 0;
        my $status=$? >> 8;
        my $locus=delete $children{$pid};
        die "Failed to download locus '$locus'. Previous database was not modified.\n" if $status != 0;
        $done++;
        print "  loci: $done/$total\r";
    }
    print "\n";
}

sub fetch_url {
    my ($url,$out)=@_;
    my $tmp="$out.part.$$";
    unlink $tmp if -e $tmp;

    if ($url =~ m{^https?://rest\.pubmlst\.org/}i) {
        # Follow the original Pneumo-Typer design: all protected PubMLST data are
        # fetched through the OAuth helper.  access_token/session_token are GLOBAL
        # in one auth directory and are not tied to a species/scheme database.
        my @cmd = (
            $^X, "$FindBin::Bin/RESTful_API_download_pubmlst.pl",
            '--url', $url,
            '--output', $tmp,
            '--auth-dir', $opt{auth_dir},
        );
        my $rc = system(@cmd);
        die "Authenticated PubMLST download failed ($url), exit=".($rc >> 8)."\n" if $rc != 0;
    }
    elsif ($url =~ m{^https?://bigsdb\.pasteur\.fr/}i) {
        die <<"PASTEUR";
The selected scheme is hosted by BIGSdb-Pasteur rather than PubMLST:
  $url
The OAuth consumer key/access_token supplied by the original Pneumo-Typer package is a
PubMLST client credential and must not be silently reused for a different BIGSdb site.
BIGSdb-Pasteur currently requires its own authenticated API access.  Select a PubMLST
scheme, or configure Pasteur credentials separately before updating this scheme.
PASTEUR
    }
    else {
        die "Unsupported MLST database host in URL: $url\n";
    }

    die "Downloaded an empty file from $url\n" unless -s $tmp;
    rename($tmp,$out) or die "Cannot rename $tmp to $out: $!\n";
}

sub slurp { my ($f)=@_; open my $fh,'<',$f or die $!; local $/; my $x=<$fh>; close $fh; return $x; }

sub validate_profiles {
    my ($f,$loci)=@_;
    open my $fh,'<',$f or die $!;
    my $h=<$fh>; close $fh;
    defined $h or die "Profiles download is empty.\n";
    chomp $h; $h =~ s/\r$//; my %x=map {$_=>1} split /\t/,$h,-1;
    my @missing=grep {!$x{$_}} @$loci;
    die "Profiles table is missing loci: ".join(', ',@missing)."\n" if @missing;
}

sub usage {
return <<'USAGE';
Usage:
  perl update_mlst_database.pl --species "Streptococcus pneumoniae"
  perl update_mlst_database.pl --species "Acinetobacter baumannii" --scheme abaumannii_2
  perl update_mlst_database.pl --list-species
  perl update_mlst_database.pl --species "Acinetobacter baumannii" --list-schemes

Options:
  --species LATIN_NAME       Species name exactly as listed in data/mlst_schemes_all.tab
  --scheme SCHEME            Scheme key; required non-interactively when one species has >1 scheme
  --scheme-description TEXT  Alternative exact scheme-description selector
  --db-root DIR              Database directory (default: ./database)
  --threads N                Parallel downloads; capped at 4 (default: 4)
  --auth-dir DIR             ONE global PubMLST OAuth token directory (default: ./auth)
  --list-species             List available species entries
  --list-schemes             List schemes for --species
  --force                    Reserved for compatibility; updates already replace only after success
USAGE
}
