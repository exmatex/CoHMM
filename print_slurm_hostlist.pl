#!/usr/bin/perl -w
use strict;
use Slurm;

die "SLURM_JOB_NODELIST is empty\n" unless defined($ENV{'SLURM_JOB_NODELIST'});
my $hl = Slurm::Hostlist::create($ENV{'SLURM_JOB_NODELIST'});
while(my $host = $hl->shift()) {
  print $host, "\n";
}
