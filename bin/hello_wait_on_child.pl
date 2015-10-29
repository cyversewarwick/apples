use strict;
use warnings;
 
use POSIX;# ":sys_wait_h";
use Time::HiRes qw(sleep);

use MooseX::Declare;
 
my $pid = fork();
my $timeout = 3;

die "Could not fork\n" if not defined $pid;

my $start_time = time();

if (not $pid) {
    print "In child\n";
    sleep 4;
    print "child finished.\n";
    exit 3;
} else {
 
    print "In parent of $pid";
    while (1) {
        my $res = waitpid($pid, WNOHANG);
        print "Res: $res\n";
        sleep(0.1);
     
        if ($res == -1) {
            print "Some error occurred ". ($? >> 8) ."\n";
            exit();
        }
        if ($res) {
            print "Child $res ended with ". ($? >> 8) ."\n";
            last;
        }

        if (time() - $start_time > $timeout) {
            print "Timeout, killing pid: $pid.\n";
            kill (SIGKILL, $pid); # SIGKILL = 15
            last;
        }
    }

}

print "Parent carries on.\n";

print "about to wait()\n";
print wait() ."\n";
print "wait() done\n";