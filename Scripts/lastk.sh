#!/bin/bash

for var in "$@"
do
  echo $var
  tail -n 5000 $var | grep "session" | tail -n 1
done
