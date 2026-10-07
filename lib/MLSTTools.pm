package MLSTTools;
use strict;
use warnings;
use Exporter 'import';
use File::Basename qw(basename);
use File::Path qw(make_path);
use POSIX qw(strftime);

our @EXPORT_OK = qw(
  load_scheme_table select_scheme list_species list_schemes
  detect_sequence_format normalize_genome safe_name which_cmd
  read_loci read_allele_index read_profiles current_time
  shell_quote
);

sub current_time { return strftime('%Y-%m-%d %H:%M:%S', localtime); }

sub load_scheme_table {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot open scheme table $file: $!\n";
    my $header = <$fh>;
    defined $header or die "Empty scheme table: $file\n";
    chomp $header;
    my @h = split /\t/, $header, -1;
    my @rows;
    while (<$fh>) {
        chomp; s/\r$//;
        next if /^\s*$/;
        my @v = split /\t/, $_, -1;
        my %r;
        @r{@h} = @v;
        push @rows, \%r;
    }
    close $fh;
    return \@rows;
}

sub list_species {
    my ($rows) = @_;
    my %x;
    for my $r (@$rows) { $x{$r->{species}}++ if defined $r->{species} && $r->{species} ne ''; }
    return [ map { [$_, $x{$_}] } sort { lc($a) cmp lc($b) } keys %x ];
}

sub list_schemes {
    my ($rows, $species) = @_;
    my @m = grep { lc($_->{species} // '') eq lc($species // '') } @$rows;
    return \@m;
}

sub select_scheme {
    my (%arg) = @_;
    my $rows    = $arg{rows} or die "rows required\n";
    my $species = $arg{species};
    my $scheme  = $arg{scheme};
    my $desc    = $arg{description};
    my @m;

    if (defined $scheme && $scheme ne '') {
        @m = grep { lc($_->{scheme} // '') eq lc($scheme) } @$rows;
        if (defined $species && $species ne '') {
            @m = grep { lc($_->{species} // '') eq lc($species) } @m;
        }
        die "Scheme '$scheme' not found" . (defined $species ? " for species '$species'" : '') . ".\n" unless @m;
    } else {
        die "--species is required unless --scheme is supplied.\n" unless defined $species && $species ne '';
        @m = grep { lc($_->{species} // '') eq lc($species) } @$rows;
        if (!@m) {
            my @suggest = grep { index(lc($_->{species} // ''), lc($species)) >= 0 || index(lc($species), lc($_->{species} // '')) >= 0 } @$rows;
            my %s; @suggest = grep { !$s{$_->{species}}++ } @suggest;
            my $msg = "Species '$species' was not found in mlst_schemes_all.tab.";
            $msg .= " Similar entries: " . join('; ', map { $_->{species} } @suggest[0 .. ($#suggest < 9 ? $#suggest : 9)]) if @suggest;
            die "$msg\n";
        }
    }

    if (defined $desc && $desc ne '') {
        my @d = grep { lc($_->{scheme_description} // '') eq lc($desc) } @m;
        die "No scheme with description '$desc' was found for '$species'.\n" unless @d;
        @m = @d;
    }

    return $m[0] if @m == 1;

    if (-t STDIN && $arg{interactive} // 1) {
        print STDERR "Multiple MLST schemes are available:\n";
        for my $i (0..$#m) {
            printf STDERR "  %d) %s | %s | %s\n", $i+1, $m[$i]{scheme}, $m[$i]{scheme_description}, $m[$i]{database};
        }
        print STDERR "Select scheme number [1-" . scalar(@m) . "]: ";
        my $ans = <STDIN>;
        chomp $ans if defined $ans;
        if (defined $ans && $ans =~ /^\d+$/ && $ans >= 1 && $ans <= @m) { return $m[$ans-1]; }
        die "Invalid scheme selection.\n";
    }

    die "Multiple schemes match. Please specify --scheme. Choices: " .
        join(', ', map { $_->{scheme} . ' [' . $_->{scheme_description} . ']' } @m) . "\n";
}

sub detect_sequence_format {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot open input $file: $!\n";
    while (<$fh>) {
        next if /^\s*$/;
        close $fh;
        return 'fasta'   if /^\s*>/;
        return 'genbank' if /^\s*LOCUS\b/;
        last;
    }
    close $fh;
    die "Cannot recognize sequence format for '$file'. Expected FASTA or GenBank/GBK.\n";
}

sub safe_name {
    my ($s) = @_;
    $s = basename($s);
    $s =~ s/[^A-Za-z0-9_.-]+/_/g;
    $s =~ s/^\.+/_/;
    return $s || 'sample';
}

sub _write_fasta_record {
    my ($fh, $id, $seq) = @_;
    $seq =~ s/[^A-Za-z]//g;
    $seq = uc($seq);
    return if $seq eq '';
    $id =~ s/\s+/_/g;
    $id =~ s/[^A-Za-z0-9_.:-]+/_/g;
    print $fh ">$id\n";
    for (my $i=0; $i<length($seq); $i+=60) { print $fh substr($seq,$i,60), "\n"; }
}

sub normalize_genome {
    my ($input, $output) = @_;
    my $fmt = detect_sequence_format($input);
    open my $out, '>', $output or die "Cannot write $output: $!\n";

    if ($fmt eq 'fasta') {
        open my $in, '<', $input or die "Cannot open $input: $!\n";
        my ($id, $seq, $n) = ('','',0);
        while (<$in>) {
            chomp;
            if (/^>(.*)/) {
                _write_fasta_record($out, $id, $seq) if $id ne '';
                $n++;
                ($id = $1) =~ s/^\s+|\s+$//g;
                $id = (split /\s+/, $id)[0] || "contig_$n";
                $seq = '';
            } else {
                s/\s+//g;
                $seq .= $_;
            }
        }
        _write_fasta_record($out, $id, $seq) if $id ne '';
        close $in;
    } else {
        open my $in, '<', $input or die "Cannot open $input: $!\n";
        my ($id, $seq, $in_origin, $record) = ('','',0,0);
        while (<$in>) {
            if (/^LOCUS\s+(\S+)/) { $id = $1; $record++; }
            elsif (/^ACCESSION\s+(\S+)/ && $id eq '') { $id = $1; }
            elsif (/^ORIGIN\b/) { $in_origin = 1; $seq = ''; }
            elsif ($in_origin && m{^//}) {
                $id ||= 'record_' . ($record || 1);
                _write_fasta_record($out, $id, $seq);
                ($id,$seq,$in_origin) = ('','',0);
            } elsif ($in_origin) {
                my $x = $_;
                $x =~ s/[^A-Za-z]//g;
                $seq .= $x;
            }
        }
        if ($in_origin && $seq ne '') {
            $id ||= 'record_' . ($record || 1);
            _write_fasta_record($out, $id, $seq);
        }
        close $in;
    }
    close $out;
    -s $output or die "No nucleotide sequence was extracted from '$input'.\n";
    return $fmt;
}

sub which_cmd {
    my ($cmd) = @_;
    for my $d (split /:/, $ENV{PATH} // '') {
        my $p = "$d/$cmd";
        return $p if -x $p && !-d $p;
    }
    return;
}

sub read_loci {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot open $file: $!\n";
    my @x;
    while (<$fh>) { chomp; s/\r$//; next if /^\s*$/; push @x, $_; }
    close $fh;
    return \@x;
}

sub read_allele_index {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot open $file: $!\n";
    my $h = <$fh>;
    my %idx;
    while (<$fh>) {
        chomp; s/\r$//;
        my ($blast_id,$locus,$allele,$original) = split /\t/, $_, 4;
        $idx{$blast_id} = [$locus,$allele,$original];
    }
    close $fh;
    return \%idx;
}

sub read_profiles {
    my ($file, $loci) = @_;
    open my $fh, '<', $file or die "Cannot open $file: $!\n";
    my $line = <$fh>;
    defined $line or die "Empty profiles file: $file\n";
    $line =~ s/^\x{FEFF}//;
    chomp $line; $line =~ s/\r$//;
    my @h = split /\t/, $line, -1;
    my %pos; for my $i (0..$#h) { $pos{$h[$i]}=$i; }
    my @missing = grep { !exists $pos{$_} } @$loci;
    die "Profiles file does not contain loci: " . join(', ', @missing) . "\n" if @missing;
    my $pk_i = 0;
    my %profiles;
    while (<$fh>) {
        chomp; s/\r$//;
        next if /^\s*$/;
        my @v = split /\t/, $_, -1;
        my @a = map { defined $v[$pos{$_}] ? $v[$pos{$_}] : '' } @$loci;
        my $key = join("\t", @a);
        $profiles{$key} = $v[$pk_i] if !exists $profiles{$key};
    }
    close $fh;
    return (\%profiles, $h[$pk_i]);
}

sub shell_quote {
    my ($s) = @_;
    $s =~ s/'/'"'"'/g;
    return "'$s'";
}

1;
