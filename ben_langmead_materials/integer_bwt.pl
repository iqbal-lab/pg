#!/usr/bin/perl -w
use strict;

my $word = shift;

$word =~ s/A/1/g;
$word =~ s/C/2/g;
$word =~ s/G/3/g;
$word =~ s/T/4/g;

sub rotations
{
    my ($str) = @_;
    my $len = length($str);
    $str = $str.$str;
    my @arr=();
    for ( 0 .. $len-1 ) 
    { 
	push @arr, substr( $str, $_, $len );       
    } 


#    print join("\n", @arr);
#    print "\n";
    return @arr;
}

sub bwm
{
    my ($str) = @_;

    my @arr= rotations($str);
    #my @arr2 = sort(@arr);
    #print join("\n",  @arr2);
    
    my @arr2 =  sort(@arr);
    return @arr2;
}

#print join("\n", bwm('abca$'));

sub bwt_via_bwm
{
    my ($str) = @_;
    my @bwm = bwm($str);
    my $i;
    my $bwt="";
    for ($i=0; $i<scalar(@bwm); $i++)
    {
	$bwt .= substr($bwm[$i],-1);
    }
    return $bwt;
}


print "bwm of $word"."\$"." is \n";
print join("\n", bwm($word.'$'));
print "\n";
#print "bwt of banana\$ is \n";
#print bwt_via_bwm('banana$');
#print "\n";


sub suffixArray
{
    my ($str) = @_;

    my @sa_tups = ();
    my $i;
    for ($i=0; $i<length($str); $i++)
    {
	my @arr = (substr($str, $i), $i);
	push @sa_tups, \@arr;
    }
}
