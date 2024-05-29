"""
This defines functions and classes to use in the annotation pipeline.
"""

import os
import subprocess
from pathlib import Path
from configparser import ConfigParser
from abc import ABC, abstractmethod
from dataclasses import dataclass
from simplify_snpeff import simplify_snpeff_default, simplify_snpeff_with_custom_annotation


# Upload a config file with the path to the SNPEff, SIFT4G and SNPEff config files
config_file_path = os.path.join(os.getcwd(), 'config.ini')
if not Path(config_file_path).exists():
    config_file_path = os.path.join(os.getcwd(), 'annotate_vcf', 'config.ini')

try:
    # Read the config file
    config = ConfigParser()
    config.read(config_file_path)

    # Rest of your code to handle the config file

except FileNotFoundError:
    print(f"Config file '{config_file_path}' not found.")


# Classes definitions:
@dataclass
class Annotation(ABC):
    """
    Abstract superclass for annotation classes.
    """
    input_file: str
    database: str
    output_folder: str

    def __post_init__(self):
        pass

    @abstractmethod
    def annotate(self):
        """
        Abstract method to perform the annotation.
        """
        pass


# SNPEff annotation subclass
@dataclass
class SNPEFFAnnotation(Annotation):
    """
    Class to annotate a vcf file with SNPEff.
    """

    # Read the config file in __init__
    snpeff: str = config.get('SNPEFF', 'snpeff')
    snpeff_config: str = config.get('SNPEFF', 'snpeff_config')

    def __post_init__(self):
        # Check if SNPEff and SNPEff config files exist
        if not Path(self.snpeff).exists():
            raise FileNotFoundError(f"SNPEff not found at {self.snpeff}")
        if not Path(self.snpeff_config).exists():
            raise FileNotFoundError(f"SNPEff config not found at {self.snpeff_config}")

        # Define input folder and prefix file
        self.input_folder, filename = os.path.split(self.input_file)
        self.prefix, _ = filename.strip().split('.vcf')

        # Check if output folder exists; create it if it doesn't
        if not Path(self.output_folder).exists():
            Path(self.output_folder).mkdir(parents=True)

        self.output_file = f"{self.output_folder}/{self.prefix}_ann.vcf"
        self.csv_stats_file = f"{self.output_folder}/{self.prefix}_ann.csv"

        self.output_file_simplified = f"{self.output_folder}/{self.prefix}_ann_simplified.vcf"

    # Annotate a vcf using database
    def annotate(self):
        """
        Method to annotate the VCF file using SNPEff.
        """
        command = (
            f"java -jar {self.snpeff} ann -c {self.snpeff_config} -classic "
            f"-noStats -csvStats {self.csv_stats_file} -v {self.database} "
            f"{self.input_file} > {self.output_file}"
        )
        # Execute the command here
        subprocess.run(command, shell=True, check=False)

    # Annotate a vcf using only custom intervals
    def annotate_with_intervals(self, interval_file: str):
        """
        Method to annotate the VCF file using SNPEff
        using a custom interval file in BED format.
        """
        command = (
            f"java -jar {self.snpeff} ann -c {self.snpeff_config} -classic "
            f"-noStats -csvStats {self.csv_stats_file} -interval {interval_file} "
            f"-v {self.database} {self.input_file} > {self.output_file}"
        )
        # Execute the command here
        subprocess.run(command, shell=True, check=False)

    # @staticmethod
    def simplify_snpeff_annotations(self, custom_annotation: str = None):
        """
        Method to simplify the SNPEff annotations.
        """
        self.output_file_simplified = f"{self.output_folder}/{self.prefix}_ann_simplified.vcf"

        # Choose the simplification method based on the custom annotation
        if custom_annotation is not None:
            simplify_snpeff_with_custom_annotation(
                file=self.output_file, custom_annotation=custom_annotation
            )

        else:
            simplify_snpeff_default(file=self.output_file)


# SIFT4G annotation subclass
@dataclass
class SIFT4GAnnotation(Annotation):
    """
    Class to annotate a vcf file with SIFT4G.
    """
    
    # Read the config file in __init__
    sift4g: str = config.get('SIFT4G', 'sift4g')

    def __post_init__(self):
        if not Path(self.sift4g).exists():
            raise FileNotFoundError(f"SIFT4G not found at {self.sift4g}")

        # Check if output folder exists; create it if it doesn't
        if not Path(self.output_folder).exists():
            Path(self.output_folder).mkdir(parents=True)

    def annotate(self):
        """
        Method to annotate the VCF file using SIFT4G.
        """
        command = (
            f"java -jar {self.sift4g} -c -i {self.input_file} -d {self.database} -r {self.output_folder}"
        )
        # Execute the command here
        subprocess.run(command, shell=True, check=False)
