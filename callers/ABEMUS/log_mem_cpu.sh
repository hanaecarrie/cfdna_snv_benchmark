#!/bin/bash

export logfile=$1
echo $logfile

touch $logfile
echo -e "date_time\tavailableRAM\tnCPU\tusageRAM\tloadaverageCPU" > $logfile


while true; do

# available RAM memory in Gb
export availRAM=$(free -h | grep Mem: | awk '{print ($2)}')
# number of CPUs
export nCPU=$(grep -c ^processor /proc/cpuinfo)
# usage RAM memory in Gb
export useRAM=$(free -h | grep Mem: | awk '{print ($3)}')
# load average
export loadavgCPU=$(uptime | awk '{print $12}' | cut -d "," -f 1)   

echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t$availRAM\t$nCPU\t$useRAM\t$loadavgCPU" >> $logfile

sleep 60
done
