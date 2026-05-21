{\rtf1\ansi\ansicpg1252\cocoartf2820
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Dental Genetics SuSiE Fine-Mapping\
\
This repository contains R workflows for SuSiE-RSS fine-mapping\
of dental GWAS traits using summary statistics and LD matrices.\
\
## Traits\
- Nteeth\
- DMFS\
- Perio\
\
## Workflow\
1. Load GWAS summary statistics\
2. Match SNPs to LD matrices\
3. Run SuSiE-RSS\
4. Extract posterior inclusion probabilities (PIPs)\
5. Summarize credible sets\
\
## Software\
- R\
- susieR\
- data.table\
- dplyr\
- DT\
\
Raw GWAS summary statistics and LD matrices are not included due to file size and/or data restrictions.}