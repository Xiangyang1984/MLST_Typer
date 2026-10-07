#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename);
use Cwd qw(abs_path);
use MLSTTools qw(load_scheme_table select_scheme list_species list_schemes normalize_genome safe_name which_cmd read_loci read_allele_index read_profiles current_time);

my %opt=(
  scheme_table=>"$FindBin::Bin/data/mlst_schemes_all.tab",
  db_root=>"$FindBin::Bin/database", auth_dir=>"$FindBin::Bin/auth",
  out=>'MLST_workplace', threads=>4, update=>'F', keep_work=>'T'
);
GetOptions(
  'input|d=s'=>\$opt{input}, 'species|s=s'=>\$opt{species}, 'scheme=s'=>\$opt{scheme},
  'scheme-description=s'=>\$opt{description}, 'output|o=s'=>\$opt{out}, 'threads|t=i'=>\$opt{threads},
  'update=s'=>\$opt{update}, 'db-root=s'=>\$opt{db_root}, 'auth-dir=s'=>\$opt{auth_dir}, 'scheme-table=s'=>\$opt{scheme_table},
  'list-species'=>\$opt{list_species}, 'list-schemes'=>\$opt{list_schemes}, 'keep-work=s'=>\$opt{keep_work},
  'help|h'=>\$opt{help}
) or die usage();
if ($opt{help}) { print usage(); exit 0; }
my $rows=load_scheme_table($opt{scheme_table});
if ($opt{list_species}) {
  print "Species\tNumber_of_schemes\n"; print $_->[0],"\t",$_->[1],"\n" for @{list_species($rows)}; exit 0;
}
if ($opt{list_schemes}) {
  die "--species is required with --list-schemes\n" unless $opt{species};
  print "database\tspecies\tscheme_description\tscheme\tURI\n";
  my $m=list_schemes($rows,$opt{species}); print join("\t",@$_{qw(database species scheme_description scheme URI)}),"\n" for @$m; exit 0;
}
die usage() unless $opt{input};
die "Input directory '$opt{input}' does not exist.\n" unless -d $opt{input};
$opt{threads}=1 if !$opt{threads} || $opt{threads}<1;
my $sel=select_scheme(rows=>$rows,species=>$opt{species},scheme=>$opt{scheme},description=>$opt{description},interactive=>1);
my $scheme=$sel->{scheme};
my $dbdir="$opt{db_root}/$scheme";
if (uc($opt{update}) eq 'T') {
  my @cmd=($^X,"$FindBin::Bin/update_mlst_database.pl",'--scheme',$scheme,'--db-root',$opt{db_root},'--scheme-table',$opt{scheme_table},'--threads',$opt{threads},'--auth-dir',$opt{auth_dir});
  push @cmd,('--species',$sel->{species});
  system(@cmd)==0 or die "MLST database update failed; typing was not started.\n";
}
for my $f (qw(alleles.fasta allele_index.tsv profiles.tsv loci.txt scheme.info.tsv)) {
  die "Database for scheme '$scheme' is incomplete ($dbdir/$f missing). Run with --update T first.\n" unless -s "$dbdir/$f";
}
my $makeblastdb=which_cmd('makeblastdb') or die "makeblastdb was not found in PATH. Install NCBI BLAST+.\n";
my $blastn=which_cmd('blastn') or die "blastn was not found in PATH. Install NCBI BLAST+.\n";

make_path($opt{out});
my $norm="$opt{out}/normalized_genomes"; my $blastdir="$opt{out}/blast"; my $dbwork="$opt{out}/blast_db";
make_path($norm,$blastdir,$dbwork);
opendir my $dh,$opt{input} or die $!;
my @files=sort grep { !/^\./ && -f "$opt{input}/$_" } readdir $dh;
closedir $dh;
die "No files were found in input directory '$opt{input}'.\n" unless @files;

my @samples;
my %used;
print "[",current_time(),"] Normalize FASTA/GenBank inputs...\n";
for my $name (@files) {
  my $safe=safe_name($name); my $base=$safe; my $n=1; while ($used{$safe}++) {$safe=$base.'.'.$n++;}
  my $nf="$norm/$safe.fasta";
  my $fmt=normalize_genome("$opt{input}/$name",$nf);
  push @samples,{name=>$name,safe=>$safe,format=>$fmt,norm=>$nf};
  print "  $name\t$fmt\n";
}

my $loci=read_loci("$dbdir/loci.txt");
my $allele_index=read_allele_index("$dbdir/allele_index.tsv");
my ($profiles,$pk_name)=read_profiles("$dbdir/profiles.tsv",$loci);

print "[",current_time(),"] MLST exact-allele search for ",scalar(@samples)," genome(s)...\n";
parallel_blast(\@samples,$dbdir,$blastdir,$dbwork,$makeblastdb,$blastn,$opt{threads});

open my $stout,'>',"$opt{out}/ST_out.txt" or die $!;
print $stout "Strain_name\tST\tProfile\n";
open my $detail,'>',"$opt{out}/MLST_detail.tsv" or die $!;
print $detail join("\t",'Strain_name','Input_format','Species','Scheme','Scheme_description','ST','Status','Profile',@$loci),"\n";
open my $long,'>',"$opt{out}/MLST_alleles.tsv" or die $!;
print $long "Strain_name\tLocus\tAllele\tCall\n";

for my $s (@samples) {
  my $calls=parse_blast("$blastdir/$s->{safe}.blastout",$allele_index,$loci);
  my (@alleles,@profile_parts); my $complete=1; my $ambig=0;
  for my $l (@$loci) {
    my $v=$calls->{$l};
    my ($allele,$call);
    if (!$v || !@$v) { ($allele,$call)=('#','missing'); $complete=0; }
    else {
      my %u=map {$_=>1} @$v; my @u=sort keys %u;
      if (@u==1) { ($allele,$call)=($u[0],'exact'); }
      else { ($allele,$call)=('?','ambiguous:'.join(',',@u)); $complete=0; $ambig=1; }
    }
    push @alleles,$allele; push @profile_parts, $l.'_'.$allele;
    print $long join("\t",$s->{name},$l,$allele,$call),"\n";
  }
  my $key=join("\t",@alleles); my ($st,$status);
  if ($complete && exists $profiles->{$key}) { $st=$profiles->{$key}; $status='Known'; }
  elsif ($complete) { $st='Novel'; $status='Novel_profile'; }
  else { $st='Unknown'; $status=$ambig ? 'Ambiguous_allele' : 'Missing_or_novel_allele'; }
  my $profile=join('-',@profile_parts);
  print $stout join("\t",$s->{name},$st,$profile),"\n";
  print $detail join("\t",$s->{name},$s->{format},$sel->{species},$scheme,$sel->{scheme_description},$st,$status,$profile,@alleles),"\n";
}
close $stout; close $detail; close $long;

open my $run,'>',"$opt{out}/run_info.tsv" or die $!;
print $run "key\tvalue\n";
for my $p (['run_time',current_time()],['input_directory',abs_path($opt{input})],['species',$sel->{species}],['scheme',$scheme],
 ['scheme_description',$sel->{scheme_description}],['source_database',$sel->{database}],['scheme_URI',$sel->{URI}],
 ['genome_number',scalar(@samples)],['locus_number',scalar(@$loci)],['profile_key_field',$pk_name]) { print $run "$p->[0]\t$p->[1]\n"; }
close $run;

if (uc($opt{keep_work}) eq 'F') { remove_tree($dbwork); }
print "[",current_time(),"] Done.\n";
print "  Main result: $opt{out}/ST_out.txt\n  Detailed result: $opt{out}/MLST_detail.tsv\n  Allele calls: $opt{out}/MLST_alleles.tsv\n";

sub parallel_blast {
  my ($samples,$dbdir,$blastdir,$dbwork,$makeblastdb,$blastn,$threads)=@_;
  $threads=1 if $threads<1; $threads=@$samples if $threads>@$samples;
  my @queue=@$samples; my %child; my $done=0; my $total=@queue;
  while (@queue || %child) {
    while (@queue && keys(%child)<$threads) {
      my $s=shift @queue; my $pid=fork(); die "fork failed: $!\n" unless defined $pid;
      if ($pid==0) {
        my $prefix="$dbwork/$s->{safe}";
        my $rc=system($makeblastdb,'-in',$s->{norm},'-dbtype','nucl','-out',$prefix);
        exit 2 if $rc!=0;
        my $out="$blastdir/$s->{safe}.blastout";
        $rc=system($blastn,'-task','blastn','-query',"$dbdir/alleles.fasta",'-db',$prefix,
          '-out',$out,'-outfmt','6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore',
          '-perc_identity','100','-qcov_hsp_perc','100','-dust','no','-max_target_seqs','100','-num_threads','1');
        exit($rc==0 ? 0 : 3);
      }
      $child{$pid}=$s->{name};
    }
    my $pid=wait(); last if $pid<0; my $status=$?>>8; my $name=delete $child{$pid};
    die "BLAST failed for '$name' (exit $status).\n" if $status!=0;
    $done++; print "  genomes: $done/$total\r";
  }
  print "\n";
}

sub parse_blast {
  my ($file,$idx,$loci)=@_;
  my %calls; my %seen;
  open my $fh,'<',$file or die "Cannot open $file: $!\n";
  while (<$fh>) {
    chomp; next if /^\s*$/;
    my @a=split /\t/;
    my ($qid,$pident,$len,$qlen,$qstart,$qend)=@a[0,2,3,4,5,6];
    next unless exists $idx->{$qid};
    next unless $pident >= 99.999 && $len == $qlen && (($qstart==1 && $qend==$qlen)||($qend==1 && $qstart==$qlen));
    my ($locus,$allele)=@{$idx->{$qid}}[0,1];
    next if $seen{"$locus\t$allele"}++;
    push @{$calls{$locus}},$allele;
  }
  close $fh;
  return \%calls;
}

sub usage {
return <<'USAGE';
Usage:
  perl mlst_typing.pl -d genomes -s "Streptococcus pneumoniae" -o result -t 8 --update T
  perl mlst_typing.pl -d genomes -s "Acinetobacter baumannii" --scheme abaumannii_2 -o result

Input:
  -d, --input DIR           Directory containing FASTA, GenBank/GBK, or a mixture of both.
                            Format is detected from file content, not filename extension.
Database selection:
  -s, --species NAME        Species name from data/mlst_schemes_all.tab.
      --scheme KEY          Select a specific scheme when a species has multiple schemes.
      --scheme-description  Select by exact scheme description.
      --update T|F          Update only the selected MLST scheme before typing (default F).
      --list-species        List available species entries.
      --list-schemes        With --species, list available schemes.
Other:
  -o, --output DIR          Output directory (default MLST_workplace).
  -t, --threads N           Number of genomes processed in parallel (default 4).
      --db-root DIR         MLST database root (default ./database).
      --auth-dir DIR        ONE shared PubMLST OAuth token directory (default ./auth).
      --keep-work T|F       Keep per-genome BLAST databases (default T).

Outputs:
  ST_out.txt                Compatible simple table: Strain_name / ST / Profile
  MLST_detail.tsv           Wide result table with species, scheme, status and per-locus calls
  MLST_alleles.tsv          Long-format per-locus allele calls
  run_info.tsv              Reproducibility metadata
USAGE
}
