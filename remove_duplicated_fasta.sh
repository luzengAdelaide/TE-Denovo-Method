#!/bin/bash

# Remove duplicated fasta sequences from Repbase library
awk '/^>/{f=!d[$1];d[$1]=1}f' in.fa > out.fa