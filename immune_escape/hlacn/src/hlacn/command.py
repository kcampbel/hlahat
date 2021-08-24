from importlib.resources import files
from hlacn import bin
from commonLib.lib.fileio import package_file_path

def bcftools_cmd(nBam, tBam, posBed, genome_fa, posTsv):
    posFormat = "'%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%AD]\n'" 
    cmd = ["bcftools", "mpileup",
            "-d", "5000", "-xBQ0", "--ignore-RG",
            "-a", "FORMAT/AD",
            "-R", posBed,
            "-f", genome_fa,
            nBam, tBam, "|",
            "bcftools", "norm", "-m-", "|",
            "bcftools", "query", "-uHf", posFormat, ">",
            posTsv]
    cmd = ' '.join(cmd)
    return(cmd)

def fitSequenza_cmd(specimen_id:str, varsFile:str, outputDir:str, altMode:str, sequenzaModelRData:str, sequenzaTools:str):
    cmd = [
        package_file_path(bin, 'fitSequenzaModel_hla.R'),
        '-i', specimen_id,
        '-v', varsFile,
        '-o', outputDir,
        '-a', altMode,
        '-r', sequenzaModelRData,
        '-t', sequenzaTools
    ]
    return(cmd)
