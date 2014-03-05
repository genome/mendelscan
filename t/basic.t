use strict;
use warnings;

use Cwd qw();
use File::Basename qw(dirname);
use File::Spec qw();
use File::Temp qw();

use Test::More tests => 3;

my %input;
subtest 'validate inputs' => sub {
    plan tests => 5;

    my $test_dir = File::Spec->join('', 'gsc', 'var', 'cache', 'testsuite', 'data', 'mendelscan', 'v1.1');
    ok(-d $test_dir, sprintf('found input directory: %s', $test_dir))
        or return;

    $input{gene} = File::Spec->join($test_dir, 'input-gene.txt');
    ok(-s $input{gene}, sprintf('found input gene file: %s', $input{gene}));

    $input{ped} = File::Spec->join($test_dir, 'input.ped');
    ok(-s $input{ped}, sprintf('found input ped file: %s', $input{ped}));

    $input{vcf} = File::Spec->join($test_dir, 'input.vcf');
    ok(-s $input{vcf}, sprintf('found input vcf file: %s', $input{vcf}));

    $input{vep} = File::Spec->join($test_dir, 'input.vep');
    ok(-s $input{vep}, sprintf('found input vep file: %s', $input{vep}));
};

my %output;
subtest 'setup outputs' => sub {
    plan tests => 3;

    my $output_dir = File::Temp->newdir();
    ok(-d $output_dir, sprintf('created output directory: %s', $output_dir))
        or return;

    $output{tsv} = File::Spec->join($output_dir, 'output.tsv');
    ok(!-e $output{tsv}, sprintf('output tsv does not exist: %s', $output{tsv}));

    $output{vcf} = File::Spec->join($output_dir, 'output.vcf');
    ok(!-e $output{vcf}, sprintf('output vcf does not exist: %s', $output{vcf}));
};

subtest 'execute MendelScan' => sub {
    plan tests => 3;

    my $parent_dir = dirname(__FILE__);
    my $top_dir = realpath($parent_dir, '..');
    my $jar = File::Spec->join($top_dir, 'dist', 'MendelScan.jar');
    ok(-s $jar, sprintf('found JAR: %s', $jar))
        or return;

    my @cmd = ('java',
        '-jar' => $jar,
        'score',
        $input{vcf},
        '--gene-file'   => $input{gene},
        '--ped-file'    => $input{ped},
        '--vep-file'    => $input{vep},
        '--output-file' => $output{tsv},
        '--output-vcf'  => $output{vcf},
    );
    system(@cmd);
    diag(join(' ', @cmd));
    my $exit = $? >> 8;
    my $signal = $? & 127;
    is($exit, 0, 'command exited zero');
    is($signal, 0, 'command did not signal');
};

sub realpath {
    return Cwd::realpath(File::Spec->join(@_));
}
