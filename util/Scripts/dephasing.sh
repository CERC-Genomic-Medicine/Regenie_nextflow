#!/bin/bash
while getopts f: flag
do
    case "${flag}" in
        f) file=${OPTARG};;
  esac
done
if grep --quiet '.|.' $file ; then sed -i '/^##/! s/|/\//g'  $file ; fi
