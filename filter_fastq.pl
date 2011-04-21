#!/usr/bin/perl
# (c) 2009 Magnus Bjursell

# standard module usage
use strict;

# custom module usage
#use lib '/home/magnus.bjursell/script/modules';
#use myFileIO qw(:basic);
#use myMathStat qw(max min round);

my $dataHR = {};
my $infoHR = {};


### =============== usage ===============###
my $progName; if ( $0 =~ /([^\/]+)$/ ) { $progName = $1; }
my $usage = qq(
Usage: $progName suffix outdir fastq_file_1 fastq_file_2

Options: suffix		Output filename suffix (fastq files will be named \'filebase.[dir].suffix.fastq\')
         outdir		Output directory (must exist)
         fastq files    Input fastq files, a pair with forward/reverse reads
\n);

my $suffix = shift; die "Please give valid suffix\n\n$usage" unless length($suffix) > 0;
my $outDir = shift; $outDir =~ s/\/$//; die "Cannot find output dir $outDir\n\n$usage" unless -d $outDir;

### ============= check options ===============###

die "\nPlease supply at least 2 input files (I found $ARGV[0] and $ARGV[1])\n\n" . $usage unless -f $ARGV[0] and -f $ARGV[1];
my $minTrmRdLen = 40;
my $maxNoN      = 4;
my $maxNoQ10low = 5;
my $maxNoQ20low = 10;


### =========== main program start ============###

{
  local $| = 1;
  print STDERR "Starting analysis...\n";

  my $inFhHR = {}; my $outFhHR = {};
  if ( $ARGV[0] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/ ) {
    $dataHR->{'inFile'}{'name'} = $1;
    $dataHR->{'inFile'}{'lane'} = $2;
    $dataHR->{'inFile'}{'dir0'} = $3;
    $inFhHR->{'0'}  = myOpen($ARGV[0]);
    my $outFn = "$outDir/" . join(".", @{$dataHR->{'inFile'}}{'name', 'lane', 'dir0'}, $suffix, "fastq");
#    die "Output file exists\n" if -f $outFn;
    $outFhHR->{'0'} = myOpenRW($outFn);
    printf STDERR ("File 1: %s\n   - Name: %s; Lane: %s; Direction: %s;\n\n", $ARGV[0], @{$dataHR->{'inFile'}}{'name', 'lane', 'dir0'});
  }

  if ( $ARGV[1] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/ ) {
    die "Multiple input samples/runs/lanes not supported\n" unless $dataHR->{'inFile'}{'name'} eq $1;
    die "Multiple input samples/runs/lanes not supported\n" unless $dataHR->{'inFile'}{'lane'} eq $2;
    $dataHR->{'inFile'}{'dir1'} = $3;
    $inFhHR->{'1'}  = myOpen($ARGV[1]);
    my $outFn = "$outDir/" . join(".", @{$dataHR->{'inFile'}}{'name', 'lane', 'dir1'}, $suffix, "fastq");
#    die "Output file exists\n" if -f $outFn;
    $outFhHR->{'1'} = myOpenRW($outFn);
    printf STDERR ("File 2: %s\n   - Name: %s; Lane: %s; Direction: %s;\n\n", $ARGV[1], @{$dataHR->{'inFile'}}{'name', 'lane', 'dir1'});
  }
  $outFhHR->{'S'} = myOpenRW("$outDir/" . join(".", @{$dataHR->{'inFile'}}{'name', 'lane'}, "S", $suffix, "fastq"));

  while ( not eof($inFhHR->{'0'}) and not eof($inFhHR->{'1'}) ) {
    my $readsHR = {};

    print STDERR "." if $. % 100000 <= 3;
    $readsHR->{'0'} = readOneSeqEntery($inFhHR->{'0'});
    $readsHR->{'1'} = readOneSeqEntery($inFhHR->{'1'});

#next if $. < 100000;
#last if $. > 400000;

    die "Error: Fastq file read pairs not matching\n" unless $readsHR->{'0'}{'nameNoDir'} eq $readsHR->{'1'}{'nameNoDir'};
    die "Error: Cannot read fastq sequence ($readsHR->{'0'}{'name'}), file $ARGV[0] probably corrupted\n" unless $readsHR->{'0'}{'name'} and $readsHR->{'0'}{'seq'} and $readsHR->{'0'}{'qual'};
    die "Error: Cannot read fastq sequence ($readsHR->{'1'}{'name'}), file $ARGV[1] probably corrupted\n" unless $readsHR->{'1'}{'name'} and $readsHR->{'1'}{'seq'} and $readsHR->{'1'}{'qual'};

    analyseOneSequence($readsHR->{'0'}, $infoHR);
    analyseOneSequence($readsHR->{'1'}, $infoHR);

    if ( $readsHR->{'0'}{'pass'} and $readsHR->{'1'}{'pass'} ) {
      my $seqDir = 0; my $outDir = 0;
      printf {$outFhHR->{$outDir}} ("%s\n%s\n+\n%s\n", @{$readsHR->{$seqDir}}{'name', 'seq', 'qual'});
      $infoHR->{'no_reads'}{'passedDir'}{$outDir}++; $infoHR->{'no_reads'}{'passedDir'}{'total'}++;
      $infoHR->{'no_bases'}{'passedDir'}{$outDir} += $readsHR->{$seqDir}{'len'}; $infoHR->{'no_bases'}{'passedDir'}{'total'} += $readsHR->{$seqDir}{'len'};

      my $seqDir = 1; my $outDir = 1;
      printf {$outFhHR->{$outDir}} ("%s\n%s\n+\n%s\n", @{$readsHR->{$seqDir}}{'name', 'seq', 'qual'});
      $infoHR->{'no_reads'}{'passedDir'}{$outDir}++; $infoHR->{'no_reads'}{'passedDir'}{'total'}++;
      $infoHR->{'no_bases'}{'passedDir'}{$outDir} += $readsHR->{$seqDir}{'len'}; $infoHR->{'no_bases'}{'passedDir'}{'total'} += $readsHR->{$seqDir}{'len'};

    } elsif ( $readsHR->{'0'}{'pass'} ) {
      my $seqDir = 0; my $outDir = "S";
      printf {$outFhHR->{$outDir}} ("%s\n%s\n+\n%s\n", @{$readsHR->{$seqDir}}{'name', 'seq', 'qual'});
      $infoHR->{'no_reads'}{'passedDir'}{$outDir}++; $infoHR->{'no_reads'}{'passedDir'}{'total'}++;
      $infoHR->{'no_bases'}{'passedDir'}{$outDir} += $readsHR->{$seqDir}{'len'}; $infoHR->{'no_bases'}{'passedDir'}{'total'} += $readsHR->{$seqDir}{'len'};

    } elsif ( $readsHR->{'1'}{'pass'} ) {
      my $seqDir = 1; my $outDir = "S";
      printf {$outFhHR->{$outDir}} ("%s\n%s\n+\n%s\n", @{$readsHR->{$seqDir}}{'name', 'seq', 'qual'});
      $infoHR->{'no_reads'}{'passedDir'}{$outDir}++; $infoHR->{'no_reads'}{'passedDir'}{'total'}++;
      $infoHR->{'no_bases'}{'passedDir'}{$outDir} += $readsHR->{$seqDir}{'len'}; $infoHR->{'no_bases'}{'passedDir'}{'total'} += $readsHR->{$seqDir}{'len'};

    }
    undef($readsHR);
  }

  foreach my $fh ( keys %{$outFhHR} ) { close($outFhHR->{$fh}); }
  foreach my $fh ( keys %{$inFhHR} )  { close($inFhHR->{$fh}); }
  print STDERR "\nAll done\n\n";
}




{ # Print information
  die "\nNo input reads\n\n" unless $infoHR->{'no_reads'}{'total'};

  print "Settings\n";
  print "Min trimmed read length\t$minTrmRdLen\n";
  print "Max number of N\t$maxNoN\n";
  print "Max number of < Q10 bases\t$maxNoQ10low\n";
  print "Max number of < Q20 bases\t$maxNoQ20low\n";
  print "\n";

  print "Summary\n";
  print "Total number of reads\t$infoHR->{'no_reads'}{'total'}\n";
  print "\n";

  printf ("Number of reads skipped because trimmed read length < %d bases\t%d\t%0.2f%%\n", $minTrmRdLen, $infoHR->{'no_reads'}{'too_short_trm_rd'}, 100 * $infoHR->{'no_reads'}{'too_short_trm_rd'} / $infoHR->{'no_reads'}{'total'});
  printf ("Number of reads skipped because of > %d 'N'\t%d\t%0.2f%%\n", $maxNoN,  $infoHR->{'no_reads'}{'too_many_N'}, 100 * $infoHR->{'no_reads'}{'too_many_N'} / $infoHR->{'no_reads'}{'total'});
  printf ("Number of reads skipped because of > %d bases under q10\t%d\t%0.2f%%\n", $maxNoQ10low, $infoHR->{'no_reads'}{'too_many_under10'}, 100 * $infoHR->{'no_reads'}{'too_many_under10'} / $infoHR->{'no_reads'}{'total'});
  printf ("Number of reads skipped because of > %d bases under q20\t%d\t%0.2f%%\n", $maxNoQ20low, $infoHR->{'no_reads'}{'too_many_under20'}, 100 * $infoHR->{'no_reads'}{'too_many_under20'} / $infoHR->{'no_reads'}{'total'});
  print "\n\n";

  die "\nNo passed reads\n\n" unless $infoHR->{'no_reads'}{'passedDir'}{'total'};

  print "Number of reads and bases passed (S: single)\n";
  print "Dir\tReads\t% of passed reads\t% of total reads\tBases\t% of passed bases\t% of total bases\n";
  foreach my $outDir ( sort { ( $a =~ /total/ ? "zzz" : $a ) cmp ( $b =~ /total/ ? "zzz" : $b ) } keys %{$infoHR->{'no_reads'}{'passedDir'}} ) {
    printf ("%s\t%d\t%0.2f%%\t%0.2f%%\t%d\t%0.2f%%\t%0.2f%%\n", $outDir, $infoHR->{'no_reads'}{'passedDir'}{$outDir}, 100 * $infoHR->{'no_reads'}{'passedDir'}{$outDir} / $infoHR->{'no_reads'}{'passedDir'}{'total'},
                   100 * $infoHR->{'no_reads'}{'passedDir'}{$outDir} / $infoHR->{'no_reads'}{'total'}, $infoHR->{'no_bases'}{'passedDir'}{$outDir}, 100 * $infoHR->{'no_bases'}{'passedDir'}{$outDir} / $infoHR->{'no_bases'}{'passedDir'}{'total'},
                   100 * $infoHR->{'no_bases'}{'passedDir'}{$outDir} / $infoHR->{'no_bases'}{'total'});
  }
  print "\n\n";

  print "Read length histogram\n";
  for my $len ( 0 .. scalar(@{$infoHR->{'length_histogram'}}) - 1 ) {
    printf ("\t%d\t%d\t%0.2f%%\n", $len, $infoHR->{'length_histogram'}[$len], 100 * $infoHR->{'length_histogram'}[$len] / $infoHR->{'no_reads'}{'passedDir'}{'total'}) if $infoHR->{'length_histogram'}[$len];
  }
  print "\n\n";

  my $smooth = 4;
  foreach my $GCpct ( sort { $a <=> $b } keys %{$infoHR->{'gc_pct_histogram'}} ) { $infoHR->{'gc_pct_histogram_avg'}{$smooth * int( ($GCpct < 100 ? $GCpct : 99.99999) / $smooth)} += $infoHR->{'gc_pct_histogram'}{$GCpct}; }
  print "Percent GC histogram\n";
  print "Percent GC\tRead count\tPercent of all reads\n";
  foreach my $GCpct ( sort { $a <=> $b } keys %{$infoHR->{'gc_pct_histogram_avg'}} ) {
    printf ("%0.2f\t%d\t%0.2f%%\n", $GCpct + ($smooth / 2), $infoHR->{'gc_pct_histogram_avg'}{$GCpct}, (100 * $infoHR->{'gc_pct_histogram_avg'}{$GCpct} / $infoHR->{'no_reads'}{'passedDir'}{'total'}) );
  }
  print "\n\n";

  print "Positionbased data\n";
  print "\tposition\tavg quality\tno 'N'\tpct 'N'\tTotal bases\n";
  for my $pos ( 1 .. scalar(@{$infoHR->{'total_bases_positionbased'}}) - 1 ) {
    printf ("\t%d\t%0.2f\t%d\t%0.2f%%\t%d\n", $pos, $infoHR->{'total_base_quality_positionbased'}[$pos] / $infoHR->{'total_bases_positionbased'}[$pos], $infoHR->{'no_N_positionbased'}[$pos],
                   100 * $infoHR->{'no_N_positionbased'}[$pos] / $infoHR->{'total_bases_positionbased'}[$pos], $infoHR->{'total_bases_positionbased'}[$pos]);
  }
  print "\n\n";

}

exit;


# SUBS

sub readOneSeqEntery {
  my $inFH = shift;
  my $readHR = {};

# Read name
  do { $readHR->{'name'} = <$inFH>; chomp($readHR->{'name'}); } until ( $readHR->{'name'} =~ /^\@/ or eof($inFH) );
  return undef if length($readHR->{'name'}) == 0;
  ($readHR->{'nameNoDir'}) = ($readHR->{'name'} =~ /([^\/]+\/)/);
# Read seq
  $readHR->{'seq'}  = <$inFH>; chomp($readHR->{'seq'});
# Check second name
  my $line = <$inFH>; chomp($line);
  die "Fastq file error: expected '+'\n" unless $line =~ /^\+/;
# Read qual
  $readHR->{'qual'} = <$inFH>; chomp($readHR->{'qual'});
  return $readHR;
}


sub analyseOneSequence {
  my $readHR = shift;
  my $infoHR = shift;

  $infoHR->{'no_reads'}{'total'}++;
  $infoHR->{'no_bases'}{'total'} += length($readHR->{'seq'});

  if ( $readHR->{'qual'} =~ /(B+)$/ ) { my $noBs = length($1); $readHR->{'seq'} = substr($readHR->{'seq'}, 0, -1 * $noBs); $readHR->{'qual'} = substr($readHR->{'qual'}, 0, -1 * $noBs); }
  $readHR->{'len'} = length($readHR->{'seq'});
  $infoHR->{'no_reads'}{'too_short_trm_rd'}++, return undef() unless $readHR->{'len'} >= $minTrmRdLen;

  my $noN = $readHR->{'seq'} =~ tr/nN//; my $noQ10low = $readHR->{'qual'} =~ tr/\@A-I//; my $noQ20low = $readHR->{'qual'} =~ tr/\@A-S//;
  $infoHR->{'no_reads'}{'too_many_N'}++, return undef() if $noN > $maxNoN;
  $infoHR->{'no_reads'}{'too_many_under10'}++, return undef() if $noQ10low > $maxNoQ10low;
  $infoHR->{'no_reads'}{'too_many_under20'}++, return undef() if $noQ20low > $maxNoQ20low;

  warn "\nNumber of bases ($readHR->{'len'}) != number of quality values (" . length($readHR->{'qual'}) . ") for $readHR->{'name'}" unless $readHR->{'len'} == length($readHR->{'qual'});
  $infoHR->{'length_histogram'}[$readHR->{'len'}]++;
  $infoHR->{'gc_pct_histogram'}{round(100 * ( $readHR->{'seq'} =~ tr/[GCgc]// ) / $readHR->{'len'}, 1)}++;

  while ( $readHR->{'seq'} =~ /[Nn]/g ) { $infoHR->{'no_N_positionbased'}[pos($readHR->{'seq'})]++; }
  while ( $readHR->{'qual'} =~ /(.)/g ) { $infoHR->{'total_base_quality_positionbased'}[pos($readHR->{'qual'})] += ord($1) - 64; $infoHR->{'total_bases_positionbased'}[pos($readHR->{'qual'})]++; }

  $readHR->{'pass'} = 1;
  return 1;
}

sub myOpen {
  my $fn = shift;
  if ( $fn eq "-" or $fn eq "STDIN" ) {
    return *STDIN;
  } else {
    open(my $fh, $fn) or die "\nCannot open file $fn for reading ($!)\n\n";
    return $fh;
  }
}

sub myOpenRW {
  my $fn = shift;
  open(my $fh, ">$fn") or die "\nCannot open file $fn for writing ($!)\n\n";
  return $fh;
}

sub min {
  my @num = @_;
  my $min = $num[0];
  foreach (@num) {$min = $_ if $_ < $min}
  return $min;
}

sub max {
  my @num = @_;
  my $max = $num[0];
  foreach (@num) {$max = $_ if $_ > $max}
  return $max;
}

sub round {
  my $num   = shift;
  my $noDec = shift; $noDec = 0 if not $noDec;
  if ( $num and $num =~ /(\-?\d*\.\d{$noDec})(\d)/ ) {
    my $saveNum = $1;
    my $detNum  = $2;
    if ( $num < 0 ) {
      if ( $detNum > 5 )  { $saveNum -= 10 ** (-1 * $noDec); }
    } else {
      if ( $detNum >= 5 ) { $saveNum += 10 ** (-1 * $noDec); }
    }
    $num = $saveNum;
  }
  return sprintf("%0." . $noDec . "f",1 * $num);
}




