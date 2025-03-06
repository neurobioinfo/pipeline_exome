#! /bin/sh
#
# $Id: splitup.sh,v 1.2 2006/05/24 23:35:50 cmhall Exp $
#
# Charles Hall <cmhall@hawaii.edu>
#
# Split a file into smaller files.  The splits are put in a directory named
# 'splits'.  The result is checked at the end by doing a diff of the original
# file and the catenation of all the splits.
#
# example:
#
#  > ./splitup.sh splitup.sh 10
#  count: 86
#  > ls -l splits
#  total 36
#  -rw-r--r--  1 cmhall users 274 2005-12-21 20:54 splitup.sh.000000000
#  -rw-r--r--  1 cmhall users 181 2005-12-21 20:54 splitup.sh.000000001
#  -rw-r--r--  1 cmhall users 539 2005-12-21 20:54 splitup.sh.000000002
#  -rw-r--r--  1 cmhall users 203 2005-12-21 20:54 splitup.sh.000000003
#  -rw-r--r--  1 cmhall users 162 2005-12-21 20:54 splitup.sh.000000004
#  -rw-r--r--  1 cmhall users 179 2005-12-21 20:54 splitup.sh.000000005
#  -rw-r--r--  1 cmhall users 229 2005-12-21 20:54 splitup.sh.000000006
#  -rw-r--r--  1 cmhall users 325 2005-12-21 20:54 splitup.sh.000000007
#  -rw-r--r--  1 cmhall users 119 2005-12-21 20:54 splitup.sh.000000008
#
# similarly, could use:
# mkdir splits
# csplit -k -z -n 9 -f splits/splitup.sh. splitup.sh 10 {*}
#
# or even better:
# split --suffix-length=9 --numeric-suffixes --lines=10 splitup.sh splitup.sh.
#
# The only downside to those is they may not behave exactly the same on all
# systems, but tail and head as used here are fairly standard.
#
# 3/23/06
# Decided to use above since it is way faster.


# check for proper number of parameters
if [ $# -lt 2 ]; then
    echo "usage: $0 <file> <step>"
    exit 1
fi

# make sure $1 is a file and it exists
if [ ! -e $1 ]; then
    echo "'$1' doesn't exist"
    exit 1
fi

# make sure $2 is a number
echo $2 | egrep "^[0-9]+$" > /dev/null || {
    echo "'step' must be a number"
    exit 1
}

# make sure $2 is not zero
if [ $2 -eq 0 ]; then
    echo "'step' must be greater than zero"
    exit 1
fi

mkdir splits || exit 1
count=`wc -l $1 | awk '{print $1}'`
echo "count: $count"

cd splits
split --suffix-length=9 --numeric-suffixes --lines=$2 ../$1 $1.
cd ..

# check the result, there should be no difference between the original
# file and the catenated splits
check=`cat splits/* | diff $1 -`
[ "$check" = "" ] || {
    echo "***********************"
    echo "something went wrong!!!"
    echo "***********************"
    exit 1
}

