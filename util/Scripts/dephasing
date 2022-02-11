#!/bin/bash

#Use : to test and unphase if necessary vcf file produced during the file generation for step 1
if grep --quiet '.|.' all.common_independent_snps.vcf ; then sed -i '/^##/! s/|/\//g'  all.common_independent_snps.vcf ; fi
