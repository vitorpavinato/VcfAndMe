"""
Script to parse config.ini file containint the path to the
SNPEff, SIFT4G and SNPEff config file.
"""

import configparser

# Instantiate the parser
config = configparser.ConfigParser()

# Parse the config.ini file fields needed
config["SNPEFF"] = {"snpeff": "/Users/tur92196/snpEff/snpEff.jar",
                    "snpeff_config": "/Users/tur92196/snpEff/snpEff.config"}

config["SIFT4G"] = {"sift4g": "/Users/tur92196/local/sift4g/SIFT4G_Annotator.jar"}

with open("annotate_vcf/config.ini", 'w', encoding="utf-8") as configfile:
    config.write(configfile)
