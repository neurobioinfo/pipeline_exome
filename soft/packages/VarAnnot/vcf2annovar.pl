#!/usr/bin/perl

use File::Basename;
use File::Path;
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use Cwd 'abs_path';
use Data::Dumper;
use Vcf;

# Command line parameters
my (@VERB, $HELP, $MAN);
my ($INFILE, $ID, $BUILDVER, $VARIANT_CALLER, $HEADER, $EXCLUDE, $INCLUDE, $RESUMEFOLDER, $MITO);

# Global variables
my ($VERBOSE1, $VERBOSE2) = ("", "");
my $ANNOVARPATH     = dirname(abs_path($0)); # (changed on tue. 8th Oct 2013) former version was : "/RQexec/dionnela/soft/packages/VarAnnot";
my $HUMANDBPATH     = "$ANNOVARPATH/humandb/"; 
my $CANCERPATH      = "$ANNOVARPATH/humandb/Cancer/";
my $CLINVARPATH     = "$ANNOVARPATH/humandb/CLINVAR/";
my $DBSNPPATH       = "$ANNOVARPATH/humandb/dbSNP/";
my $ENCODEPATH      = "$ANNOVARPATH/humandb/Encode.UCSC/";
my $FREQUENCIESPATH = "$ANNOVARPATH/humandb/Frequencies";
my $GENANNOPATH     = "$ANNOVARPATH/humandb/Genanno/";
my $GENOMETRAXPATH  = "$ANNOVARPATH/humandb/GenomeTrax/";
my $BODYMAPPATH     = "$ANNOVARPATH/humandb/IlluminaBodyMap2/";
my $MIRNAPATH       = "$ANNOVARPATH/humandb/miRNA/";
my $MITOPATH        = "$ANNOVARPATH/humandb/Mitochondria/";
my $REPEATPATH      = "$ANNOVARPATH/humandb/REPEATS/";
my $SCORESPATH      = "$ANNOVARPATH/humandb/Scores";
my $DBNSFPPATH      = "$ANNOVARPATH/humandb/Scores/dbNSFP";
my $WG_GERPPATH     = "$ANNOVARPATH/humandb/Scores/GERP";
my $GWAVAPATH       = "$ANNOVARPATH/humandb/Scores/GWAVA";
my $DBSCSNVPATH     = "$ANNOVARPATH/humandb/Scores/dbSCSNV";
my $WGSDANNPATH     = "$ANNOVARPATH/humandb/Scores/DANN";
my $WGSCADDPATH     = "$ANNOVARPATH/humandb/Scores/CADD";
my $EXACPATH        = "$ANNOVARPATH/humandb/Frequencies/ExAC";
my $HRCR1PATH       = "$ANNOVARPATH/humandb/Frequencies/HRCR1";
my $CGPATH          = "$ANNOVARPATH/humandb/Frequencies/CG";
my $GMEPATH         = "$ANNOVARPATH/humandb/Frequencies/GME";
my $GNOMADPATH      = "$ANNOVARPATH/humandb/Frequencies/GNOMAD";
my $WGSEIGENPATH    = "$ANNOVARPATH/humandb/Scores/WGSEIGEN";
my $WGSFATHMMPATH   = "$ANNOVARPATH/humandb/Scores/FATHMM";
my $KAVIARPATH      = "$ANNOVARPATH/humandb/Frequencies/Kaviar";
my $RVIS_EXACPATH   = "$ANNOVARPATH/humandb/Scores/RVIS_ExAC";
my $REVELPATH       = "$ANNOVARPATH/humandb/Scores/REVEL";
my $INTERVARPATH    = "$ANNOVARPATH/humandb/INTERVAR";
my $MPC_EXACPATH    = "$ANNOVARPATH/humandb/Scores/MPC_ExAC";
my $UORFPATH        = "$ANNOVARPATH/humandb/uORF";
my $UVSIGNATUREPATH = "/lustre03/project/6001220/COMMON/runs/spiegelm/bradley"; # DAN: temporary hard code until Alex gives access to humandb folder
my $ANNOINFILE;                 # path of the annovar input file
my $OUTFILE;                    # path of the output file
my $OUTFILE_ANNO;               # path of the annovar --geneanno output file
my $OUTFILE_DBSNP138;           # path of the annovar --filter snp138 output file
my $OUTFILE_DBSNP138_COMMON;    # path of the annovar --filter snp138common output file
my $OUTFILE_DBSNP138_FLAGGED;   # path of the annovar --filter snp138flagged output file
my $OUTFILE_DBSNP142;           # path of the annovar --filter snp142 output file
my $OUTFILE_DBSNP142_COMMON;    # path of the annovar --filter snp142common output file
my $OUTFILE_DBSNP142_FLAGGED;   # path of the annovar --filter snp142flagged output file
my $OUTFILE_DBSNP132;           # path of the annovar --filter snp132 output file
my $OUTFILE_1000G;              # path of the annovar --filter 1000g output file
my $OUTFILE_EAS;                # path of the annovar --filter 1000g east asian allele frequency output file
my $OUTFILE_SAS;                # path of the annovar --filter 1000g south asian allele frequency output file
my $OUTFILE_AFR;                # path of the annovar --filter 1000g african allele frequency output file
my $OUTFILE_AMR;                # path of the annovar --filter 1000g ad mixed american allele frequency output file
my $OUTFILE_EUR;                # path of the annovar --filter 1000g european allele frequency output file
my $OUTFILE_EXAC;               # path of the annovar --filter exac output file
my $OUTFILE_CG46;               # path of the annovar --filter Complete Genomics output file
my $OUTFILE_GME;                # path of the annovar --filter Great Middle East Frequencies
my $OUTFILE_GNOMAD_EXOME;       # path of the annovar --filter gnomAD exome allele frequency output file
my $OUTFILE_GNOMAD_WGS;         # path of the annovar --filter gnomAD genome allele frequency output file
my $OUTFILE_GONL;               # path of the annovar --filter gonl output file
my $OUTFILE_WELLDERLY;          # path of the annovar --filter wellderly output file
my $OUTFILE_NCI60;              # path of the annovar --filter NCI60 output file
my $OUTFILE_DBNSFP;             # path of the annovar --filter dbnsfp output file
my $OUTFILE_GERP_GT2_WG;            # path of the annovar --filter gerp++gt2 output file
my $OUTFILE_GERP_ELEM_WG;       # path of the annovar --regionanno gerp++elem output file
my $OUTFILE_ESP6500;            # path of the annovar --filter esp6500 output file
my $OUTFILE_ESP6500AA;          # path of the annovar --filter esp6500aa output file
my $OUTFILE_ESP6500EA;          # path of the annovar --filter esp6500ea output file
my $OUTFILE_CLINVAR;            # path of the annovar --filter clinvar_20170130 output file
my $OUTFILE_COSMIC;             # path of the annovar --filter cosmic70 output file
my $OUTFILE_ICGC21;             # path of the annovar --filter International Cancer Genome Consortium version 21 output file
my $OUTFILE_WGSCADD;            # path of the annovar --filter cadd output file
my $OUTFILE_WGSDANN;            # path of the annovar --filter dann output file
my $OUTFILE_DBSCSNV11;          # path of the annovar --filter dbscsnv11 output file
my $OUTFILE_HRCR1;              # path of the annovar --filter hrcr1 output file
my $OUTFILE_WGSFATHMM;          # path of the annovar --filter fathmm output file
my $OUTFILE_GWAVA;              # path of the annovar --filter gwava output file
my $OUTFILE_WGSEIGEN;      	# path of the annovar --filter eigen output file
my $OUTFILE_ExAC_ALL_nonpsych;  # path of the annovar --filter exac_nonpsych output file
my $OUTFILE_ExAC_ALL_nontcga;   # path of the annovar --filter exac_nontcga output file
my $OUTFILE_KAVIAR;             # path of the annovar --filter kaviar output file
my $OUTFILE_GTX_GWAS;           # path of the annovar --filter gtx_gwas output file
my $OUTFILE_GTX_HGMD;           # path of the annovar --filter gtx_hgmd output file
my $OUTFILE_REVEL;              # path of the annovar --filter reval output file
my $OUTFILE_INTERVAR;           # path of the annovar --filter intervar output file
my $OUTFILE_EBI3222PAK;         # path of the annovar --filter EBI3222PAK output file
my $OUTFILE_GTX_CHIPSEQ;        # path of the annovar --regionanno gtx_chipseq output file
my $OUTFILE_GTX_CPG;            # path of the annovar --regionanno gtx_cpg output file
my $OUTFILE_GTX_DISEASE;        # path of the annovar --regionanno gtx_disease output file
my $OUTFILE_GTX_DNASE;          # path to the annovar --regionanno dnase_hg19 output file
my $OUTFILE_GTX_DRUG;           # path of the annovar --regionanno gtx_drug output file
my $OUTFILE_GTX_HGMD_DISGENE;   # path of the annovar --regionanno gtx_hgmd_disgene output file
my $OUTFILE_GTX_HGMDIMPUTED;    # path of the annovar --regionanno gtx_hgmdimputed output file
my $OUTFILE_GTX_MICROSAT;       # path of the annovar --regionanno gtx_microsat output file
my $OUTFILE_GTX_MIRNA;          # path of the annovar --regionanno gtx_mirna output file
my $OUTFILE_GTX_OMIM;           # path of the annovar --regionanno gtx_omim output file
my $OUTFILE_GTX_ORPHA;          # path of the annovar --regionanno gtx_orpha output file
my $OUTFILE_GTX_PATH;           # path of the annovar --regionanno gtx_path output file
my $OUTFILE_GTX_PGMD;           # path of the annovar --regionanno gtx_pgmd output file
my $OUTFILE_GTX_PTMS;           # path of the annovar --regionanno gtx_ptms output file
my $OUTFILE_GTX_TRANSFAC_TFBS;  # path of the annovar --regionanno gtx_transfac output file
my $OUTFILE_GTX_TSS;            # path of the annovar --regionanno gtx_tss output file
my $OUTFILE_GTX_VISTA;          # path of the annovar --regionanno gtx_vista output file
my $OUTFILE_ENC_HMM;            # path of the annovar --regionanno enc_hmm output file
my $OUTFILE_ENC_TFBS;           # path of the annovar --regionanno enc_tfbs output file
my $OUTFILE_SIMPLE_REPEAT;      # path of the annovar --regionanno simple_repeat output file
my $OUTFILE_SEGMENTAL_DUPS;     # path of the annovar --regionanno segment_dups output file
my $OUTFILE_REPEAT_MASKER;      # path of the annovar --regionanno repeat_masker output file
my $OUTFILE_MIRNA;              # path of the annovar --regionanno mirna output file
my $OUTFILE_MIRNATARGET;        # path of the annovar --regionanno mirnatarget output file
my $OUTFILE_RNASEQ_CUFF;        # path of the annovar --regionanno rnaseq_cuff output file
my $OUTFILE_RNASEQ_SCRI;        # path of the annovar --regionanno rnaseq_scri output file
my $OUTFILE_RVIS_EXAC;          # path of the annovar --regionanno rvis_exac output file
my $OUTFILE_MPC_EXAC;           # path of the annovar --filter rvis_mpc output file
my $OUTFILE_UORF;               # path of the annovar --filter uorf output file
my $OUTFILE_UV_SIGNATURE;       # path of the annovar --filter uv_signature output file
my $OUTPATH;                    # path of the results output file (if not provided, will be interpreted from $OUTFILE, or even $INFILE if $OUTFILE is not provided)
my $SUBOUTPATH;                 # path of the all the annovar output files (if not provided, will be interpreted from $OUTFILE, or even $INFILE if $OUTFILE is not provided)
my $USEDBUILD;                  # build to be use during the annovar analysis (can be different than $BUILDVER, which will be the one printed in name of output files)
my @EXCLUDETAB = ();
my @INCLUDETAB = ();

&main;

sub main {
    # First process the provided arguments
    processArguments();

    # Parse input data into hash
    $VERBOSE1 and printerr("\n**********\n* Starting parsing $VARIANT_CALLER input\n**********\n\n");
    my ($rH_data, $rH_vcfVarIndex) = convertVcf();
 
    $VERBOSE1 and printerr("\n**********\n* Parsing $VARIANT_CALLER input done...\n**********\n\n");

    my $rO_vcf = Vcf->new(file => $INFILE);
    $rO_vcf->parse_header();

    # Now proceed to the sequential (-geneanno, dbsnp, 1000g...) executions of Annovar
    launchAnnovar(
        -data => $rH_data,
        -vcf  => $rO_vcf
    );

    # Build object to be printed out, specifically to the wanted output format
    $VERBOSE1 and printerr("\n**********\n* Starting printing VCF output\n**********\n\n");
    printVcf(
        -data  => $rH_data,
        -vcf   => $rO_vcf,
        -index => $rH_vcfVarIndex
    );
    $VERBOSE1 and printerr("\n**********\n* Printing VCF output done...\n**********\n\n");

    print "\n\nAnnotation of $VARIANT_CALLER $BUILDVER variants for $ID is a great success !!\nThank you for using this program, see you next time !! :)\n\n\n";
}

sub processArguments {
    my @command_line = @ARGV;       # command line argument
    GetOptions('verbose|v'=>\@VERB, 'help|h'=>\$HELP, 'man|m'=>\$MAN, 'infile|i=s'=>\$INFILE, 'outfile|o=s'=>\$OUTFILE, 'identifier|id=s'=>\$ID, 'buildver|b=s'=>\$BUILDVER, 'vcaller|vc=s'=>\$VARIANT_CALLER, 'exclude|e=s'=>\$EXCLUDE, 'include|inc=s'=>\$INCLUDE, 'rfolder|rf=s'=>\$RESUMEFOLDER, 'mito'=>\$MITO) || pod2usage();

    $HELP and pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDOUT);
    $MAN  and pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDOUT);
    @ARGV == 0 || pod2usage("\n**************************************************************\n Syntax error\n\n\n");

    if ($INCLUDE && $EXCLUDE && ($EXCLUDE !~ /^all$/i && $INCLUDE !~ /^all$/)) {
        pod2usage("\n**************************************************************\nError in argument: you cannot use both --include and -exclude arguments at the same time !!\n\n\n");
    }

    if (not $ID) {
        pod2usage("\n**************************************************************\nError in argument: the --id has to be provided !!\n\n\n");
    }
    $ID = ucfirst $ID;

    if (($VARIANT_CALLER) && ($VARIANT_CALLER !~ /^(varscan|dibayes|smallindel|smallindels|mpileup|gatk|dindel|custom)$/i)) {
        pod2usage("\n**************************************************************\nError in argument: the --vcaller must be one of these : varscan, dibayes, smllindels, mpileup, gatk, dindel, custom\n\n\n");
    }
    $VARIANT_CALLER = lc $VARIANT_CALLER if ($VARIANT_CALLER);

    if ((not $BUILDVER) || ($BUILDVER !~ /^(hg19|19|b37|v37)$/i)) {
        pod2usage("\n**************************************************************\nError in argument: the --buildver must be one of these : hg19 (or 19), b37, v37\n\n\n");
    } elsif (($BUILDVER !~ /^hg/i) && ($BUILDVER !~ /b37|v37/i)) {
        $BUILDVER = "hg" . $BUILDVER;
    }
    $BUILDVER = lc $BUILDVER;
    if ($EXCLUDE) {
        ($EXCLUDE !~ /^((all|gene|snp132|snp138|snp138common|snp138flagged|snp142|snp142common|snp142flagged|1kg|1kg_afr|1kg_amr|1kg_eur|1kg_eas|1kg_sas|exac|exac_afr|exac_amr|exac_eas|exac_fin|exac_nfe|exac_oth|exac_sas|cg|gme|gnomad_exome|gnomad_genome|dbnsfp|sift|pp2|lrt|mt|ma|fathmm|provean|vest|mcap|cadd|dann|fatmmkl|msvm|mlr|intgr|gerp|gerp\+\+gt2|gerp\+\+elem|phylop|phcons|siphy|eigen|genocanyon|interpro|gtex|esp6500|esp6500aa|esp6500ea|mirna|mirnatarget|rnaseq_cuff|rnaseq_scri|clinvar|cosmic|icgc21|nci60|gtx_chipseq|gtx_cpg|gtx_disease|gtx_dnase|gtx_drug|gtx_gwas|gtx_hgmd|gtx_hgmd_disgene|gtx_hgmdimputed|gtx_microsat|gtx_mirna|gtx_omim|gtx_orpha|gtx_path|gtx_pgmd|gtx_ptms|gtx_transfac|gtx_tss|gtx_vista|enc_hmm|enc_tfbs|simple_repeat|segment_dups|rmsk|gonl|wellderly|wgsdann|gwava|wgscadd|wgsfathmm|hrcr1|dbscsnv11|exac_nonpsych|exac_nontcga|kaviar|wgseigen|rvis_exac|revel|intervar|ebi3222pak|mpc_exac|uorf|uv_signature),?)+$/i) and pod2usage("\n**************************************************************\nError in argument: the --exclude must be one of these : gene,snp132,snp138,snp138common,snp138flagged,snp142,snp142common,snp142flagged,1kg,1kg_afr,1kg_amr,1kg_eur,1kg_eas,1kg_sas,exac,exac_afr,exac_amr,exac_eas,exac_fin,exac_nfe,exac_oth,exac_sas,cg,gme,gnomad_exome,gnomad_genome,dbnsfp,sift,pp2,lrt,mt,ma,fathmm,provean,vest,cadd,mcap,dann,fatmmkl,msvm,mlr,intgr,gerp,gerp\+\+gt2,gerp\+\+elem,phylop,phcons,siphy,eigen,genocanyon,interpro,gtex,esp6500,esp6500aa,esp6500ea,mirna,mirnatarget,rnaseq_cuff,rnaseq_scri,clinvar,cosmic,icgc21,nci60,gtx_chipseq,gtx_cpg,gtx_disease,gtx_dnase,gtx_drug,gtx_gwas,gtx_hgmd,gtx_hgmd_disgene,gtx_hgmdimputed,gtx_microsat,gtx_mirna,gtx_omim,gtx_orpha,gtx_path,gtx_pgmd,gtx_ptms,gtx_transfac,gtx_tss,gtx_vista,enc_hmm,enc_tfbs,simple_repeat,segment_dups,rmsk,gonl,wellderly,wgsdann,gwava,wgscadd,wgsfathmm,hrcr1,dbscsnv11,exac_nonpsych,exac_nontcga,kaviar,wgseigen,rvis_exac,revel,intervar,ebi3222pak,mpc_exac,uorf,uv_signature\n\n\n");
        @EXCLUDETAB = split /,/, $EXCLUDE; # / just to avoid a bug in eclipse display.../
    }

    if ($INCLUDE) {
        ($INCLUDE !~ /^((all|gene|snp132|snp138|snp138common|snp138flagged|snp142|snp142common|snp142flagged|1kg|1kg_afr|1kg_amr|1kg_eur|1kg_eas|1kg_sas|exac|exac_afr|exac_amr|exac_eas|exac_fin|exac_nfe|exac_oth|exac_sas|cg|gme|gnomad_exome|gnomad_genome|dbnsfp|sift|pp2|lrt|mt|ma|fathmm|provean|vest|cadd|dann|fatmmkl|msvm|mlr|intgr|gerp|gerp\+\+gt2|gerp\+\+elem|phylop|phcons|mcap|eigen|genocanyon|interpro|gtex|siphy|esp6500|esp6500aa|esp6500ea|mirna|mirnatarget|rnaseq_cuff|rnaseq_scri|clinvar|cosmic|icgc21|nci60|gtx_chipseq|gtx_cpg|gtx_disease|gtx_dnase|gtx_drug|gtx_gwas|gtx_hgmd|gtx_hgmd_disgene|gtx_hgmdimputed|gtx_microsat|gtx_mirna|gtx_omim|gtx_orpha|gtx_path|gtx_pgmd|gtx_ptms|gtx_transfac|gtx_tss|gtx_vista|enc_hmm|enc_tfbs|simple_repeat|segment_dups|rmsk|gonl|wellderly|wgsdann|gwava|wgscadd|wgsfathmm|hrcr1|dbscsnv11|exac_nonpsych|exac_nontcga|kaviar|wgseigen|rvis_exac|revel|intervar|ebi3222pak|mpc_exac|uorf|uv_signature),?)+$/i) and pod2usage("\n**************************************************************\nError in argument: the --include must be one of these : gene,snp132,snp138,snp138common,snp138flagged,snp142,snp142common,snp142flagged,1kg,1kg_afr,1kg_amr,1kg_eur,1kg_eas,1kg_sas,exac,exac_afr,exac_amr,exac_eas,exac_fin,exac_nfe,exac_oth,exac_sas,cg,gme,gnomad_exome,gnomad_genome,dbnsfp,sift,pp2,lrt,mt,ma,fathmm,provean,vest,cadd,dann,fatmmkl,msvm,mlr,intgr,gerp,gerp++gt2,gerp++elem,phylop,phcons,mcap,eigen,genocanyon,interpro,gtex,siphy,esp6500,esp6500aa,esp6500ea,mirna,mirnatarget,rnaseq_cuff,rnaseq_scri,clinvar,cosmic,icgc21,nci60,gtx_chipseq,gtx_cpg,gtx_disease,gtx_dnase,gtx_drug,gtx_gwas,gtx_hgmd,gtx_hgmd_disgene,gtx_hgmdimputed,gtx_microsat,gtx_mirna,gtx_omim,gtx_orpha,gtx_path,gtx_pgmd,gtx_ptms,gtx_transfac,gtx_tss,gtx_vista,enc_hmm,enc_tfbs,simple_repeat,segment_dups,rmsk,gonl,wellderly,wgsdann,gwava,wgscadd,wgsfathmm,hrcr1,dbsvsnv11,exac_nonpsych,exac_nontcga,kaviar,wgseigen,rvis_exac,revel,intervar,ebi3222pak,mpc_exac,uorf,uv_signature\n\n\n");
        @INCLUDETAB = split /,/, $INCLUDE; # / just to avoid a bug in eclipse display...
        #split dbnsfp in many databases
        if ( grep( /^dbnsfp$/, @INCLUDETAB ) ) {
          push (@INCLUDETAB, qw(sift pp2 lrt mt ma fathmm provean vest msvm mlr mcap cadd dann fathmmkl eigen genocanyon intgr gerp phylop phcons siphy interpro gtex));
          @INCLUDETAB = do { my %seen; grep { !$seen{$_}++ } @INCLUDETAB };
        }
    }

    if (scalar(@VERB) == 1) {
        $VERBOSE1 = "--verbose";
    } elsif (scalar(@VERB) == 2) {
        $VERBOSE1 = "--verbose";
        
        $VERBOSE2 = "--verbose";
    }

    if (not $INFILE) {
        pod2usage("\n**************************************************************\nError in argument: the --infile has to be provided !!\n\n\n");
    }

    my $outFileNotice = undef;
    if (not $OUTFILE) {
        my ($directory, $filename) = $INFILE =~ m/(.*\/)(.*)$/;
        (not $directory) and $directory = "./";
        $OUTPATH = $directory;
    } else {
        my ($directory, $filename) = $OUTFILE =~ m/(.*\/)(.*)$/;
        (not $directory) and $directory = "./";
        if ($filename && -d $filename) {
            $directory .= $filename . "/";
            $filename = undef;
        }
        $OUTPATH = $directory;
        (not $filename) and $OUTFILE = "";
        if (not -e $OUTPATH) {
            my $print = "path $OUTPATH does not exists...";
            my ($directory, $filename) = $INFILE =~ m/(.*\/)(.*)$/;
            (not $directory) and $directory = "./";
            $OUTPATH = $directory;
            $outFileNotice = "NOTICE : $print $OUTPATH will be used instead...\n";
        }
    }
    if (not $OUTPATH) {
        my ($directory, $filename) = $INFILE =~ m/(.*\/)(.*)$/;
        (not $directory) and $directory = "./";
        $OUTPATH = $directory;
    }
    $OUTPATH .= "/" if ($OUTPATH !~ /\/$/);

    if ($RESUMEFOLDER) {
        (!(-d $RESUMEFOLDER)) and pod2usage("\n**************************************************************\nError in argument: the provided --rfolder $RESUMEFOLDER is not a folder !!\n\n\n");
        (!(-e $RESUMEFOLDER)) and pod2usage("\n**************************************************************\nError in argument: the provided --rfolder $RESUMEFOLDER does not exists !!\n\n\n");
        $RESUMEFOLDER =~ s/\/$//;
        $SUBOUTPATH = $RESUMEFOLDER . "/";
    } else {
        my $date = GetCurrentTime();
        $SUBOUTPATH = $OUTPATH . $ID . "_ANNOVAR_analysis_files_" . $date . "/";
        mkpath $SUBOUTPATH, 002;
    }

    $ANNOINFILE                 = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_input." . $VARIANT_CALLER . "." . $BUILDVER : $SUBOUTPATH . $ID . "_ANNOVAR_input." . $BUILDVER;
    $OUTFILE_ANNO               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_geneanno_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_geneanno_output";
    $OUTFILE_DBSNP132           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp132_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp132_output";
    $OUTFILE_DBSNP138           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp138_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp138_output";
    $OUTFILE_DBSNP138_COMMON    = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp138Common_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp138Common_output";
    $OUTFILE_DBSNP138_FLAGGED   = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp138Flagged_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp138Flagged_output";
    $OUTFILE_DBSNP142           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp142_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp142_output";
    $OUTFILE_DBSNP142_COMMON    = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp142Common_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp142Common_output";
    $OUTFILE_DBSNP142_FLAGGED   = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp142Flagged_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp142Flagged_output";
    $OUTFILE_1000G              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_output";
    $OUTFILE_EAS                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_eas_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_eas_output";
    $OUTFILE_SAS                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_sas_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_sas_output";
    $OUTFILE_AFR                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_afr_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_afr_output";
    $OUTFILE_AMR                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_amr_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_amr_output";
    $OUTFILE_EUR                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_eur_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_eur_output";
    $OUTFILE_EXAC               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_exac_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_exac_output";
    $OUTFILE_CG46               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_cg46_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_cg46_output";
    $OUTFILE_GME                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gme_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gme_output";
    $OUTFILE_GNOMAD_EXOME       = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gnomad_exome_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gnomad_exome_output";
    $OUTFILE_GNOMAD_WGS         = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gnomad_genome_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gnomad_genome_output";
    $OUTFILE_COSMIC             = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_cosmic_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_cosmic_output";
    $OUTFILE_ICGC21             = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_icgc21_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_icgc21_output";
    $OUTFILE_WGSCADD            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_wgscadd_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_wgscadd_output";
    $OUTFILE_WGSDANN            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_wgsdann_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_wgsdann_output"; 
    $OUTFILE_WGSEIGEN           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_WGSEIGEN_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_wgseigen_output";
    $OUTFILE_DBSCSNV11          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbscsnv11_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbscsnv11_output";
    $OUTFILE_HRCR1              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_hrcr1_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_hrcr1_output";
    $OUTFILE_WGSFATHMM          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_fathmm_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_fathmm_output";
    $OUTFILE_GWAVA              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gwava_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gwava_output";
    $OUTFILE_ExAC_ALL_nonpsych  = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_exac_all_nonpsych_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_exac_all_nonpsych_output";
    $OUTFILE_ExAC_ALL_nontcga   = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_exac_all_nontcga_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_exac_all_nontcga_output";
    $OUTFILE_KAVIAR             = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_kaviar_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_kaviar_output";
    $OUTFILE_NCI60              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_nci60_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_nci60_output";
    $OUTFILE_GONL               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gonl_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gonl_output";
    $OUTFILE_WELLDERLY          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_wellderly_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_wellderly_output";
    $OUTFILE_DBNSFP             = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbnsfp33a_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbnsfp33a_output";
    $OUTFILE_GERP_GT2_WG        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gerp++gt2_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gerp++gt2_output";
    $OUTFILE_GERP_ELEM_WG       = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gerp++elem_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gerp++elem_output.";
    $OUTFILE_ESP6500            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_esp6500_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_esp6500_output";
    $OUTFILE_ESP6500AA          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_esp6500aa_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_esp6500aa_output";
    $OUTFILE_ESP6500EA          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_esp6500ea_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_esp6500ea_output";
    $OUTFILE_MIRNA              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_mirna_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_mirna_output";
    $OUTFILE_MIRNATARGET        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_mirnatarget_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_mirnatarget_output";
    $OUTFILE_RNASEQ_CUFF        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Cufflinks_allTissues_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Cufflinks_allTissues_output";
    $OUTFILE_RNASEQ_SCRI        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Scripture_allTissues_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Scripture_allTissues_output";
    $OUTFILE_GTX_CHIPSEQ        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_chipseq_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_chipseq_output";
    $OUTFILE_CLINVAR            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_clinvar_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_clinvar_output";
    $OUTFILE_GTX_CPG            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_cpg_islands_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_cpg_islands_output";
    $OUTFILE_GTX_DISEASE        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_disease_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_disease_output";
    $OUTFILE_GTX_DNASE          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_dnase_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_dnase_output";
    $OUTFILE_GTX_DRUG           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_drug_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_drug_output";
    $OUTFILE_GTX_GWAS           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_gwas_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_gwas_output";
    $OUTFILE_GTX_HGMD           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_output";
    $OUTFILE_GTX_HGMD_DISGENE   = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_disease_genes_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_disease_genes_output";
    $OUTFILE_GTX_HGMDIMPUTED    = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmdimputed_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmdimputed_output";
    $OUTFILE_GTX_MICROSAT       = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_microsatellites_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_microsatellites_output";
    $OUTFILE_GTX_MIRNA          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_mirna_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_mirna_output";
    $OUTFILE_GTX_OMIM           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_omim_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_omim_output";
    $OUTFILE_GTX_ORPHA           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_orpha_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_orpha_output";
    $OUTFILE_GTX_PATH           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pathway_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pathway_output";
    $OUTFILE_GTX_PGMD           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pgmd_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pgmd_output";
    $OUTFILE_GTX_PTMS           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_ptms_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_ptms_output";
    $OUTFILE_GTX_TRANSFAC_TFBS  = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_transfac_sites_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_transfac_sites_output";
    $OUTFILE_GTX_TSS            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_tss_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_tss_output";
    $OUTFILE_GTX_VISTA          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_vista_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_vista_output";
    $OUTFILE_ENC_HMM            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_EncodeUCSC_hmm_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_EncodeUCSC_hmm_output";
    $OUTFILE_ENC_TFBS           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_EncodeUCSC_tbfs_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_EncodeUCSC_tbfs_output";
    $OUTFILE_SIMPLE_REPEAT      = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_simpleRepeat_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_simpleRepeat_output";
    $OUTFILE_SEGMENTAL_DUPS     = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_segmentalDups_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_segmentalDups_output";
    $OUTFILE_REPEAT_MASKER      = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_repeatMasker_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_repeatMasker_output";
    $OUTFILE_RVIS_EXAC          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_rvis_exac_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_rvis_exac_output";
    $OUTFILE_REVEL              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_revel_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_revel_output";
    $OUTFILE_INTERVAR           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_intervar_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_intervar_output";
    $OUTFILE_EBI3222PAK         = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_ebi3222pak_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_ebi3222pak_output";
    $OUTFILE_MPC_EXAC           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_mpc_exac_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_mpc_exac_output";
    $OUTFILE_UORF               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_uorf_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_uorf_output";
    $OUTFILE_UV_SIGNATURE       = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_uvSignature_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_uvSignature_output";


    my $logFile = $SUBOUTPATH . $ID . "_ANNOVAR_unresolved.log";
    open(LOG, ">$logFile") || die "Could not open log file $logFile...\n($!)\n";

    ($outFileNotice && $VERBOSE1) && printerr($outFileNotice);
}

sub convertVcf {
    if ($RESUMEFOLDER && -e "$ANNOINFILE") {
        # resuming ANNOVAR without converting vcf
        $VERBOSE1 and printerr("\n**********\n* Resuming ANNOVAR analysis using already existing converted file $ANNOINFILE\n**********\n\n");
    } else {
        # Execute convert2annovar to prepare the annovar input file
#        my $systemCommand = "perl $ANNOVARPATH/convert2annovar.pl -format vcf4old -allallele -includeinfo -outfile $ANNOINFILE $VERBOSE2 $INFILE";
        my $systemCommand = "perl $ANNOVARPATH/convert2annovar.pl -format vcf4 -includeinfo -allsample -withfreq -outfile $ANNOINFILE $VERBOSE2 $INFILE";
        $VERBOSE1 and printerr("\n**********\n* Starting converting vcf input to annovar input file with system command:\n* <$systemCommand>\n**********\n\n");
        my $error = system ($systemCommand);
        if ($error) {
            printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
        } elsif ($VERBOSE2) {
            my $exitSignal = ($? >> 8);
            my $message = sprintf "Command exited with value %d", $exitSignal;
            printerr("\n**********\n* $message\n**********\n\n");
        }
        $VERBOSE1 and printerr("\n**********\n* Input conversion done...\n**********\n\n");
    }

    open(ANNOINFILE, "<$ANNOINFILE") || die "cannot open file $ANNOINFILE:\n$!\n";

    my %vcfVariants;
    my %input;
    my @annoinfile;

    $VERBOSE1 and printerr("\n**********\n* Preprocessing of the converted file : assigning Variant Class to variants & building of the variant hash...\n**********\n\n");
    # Parse converted file to prepare VariantClass (VC), build a hash (%vcfVariants) to track back variants from original vcf input file
    # and also remove some trailing columns to make it lighter.
    while (my $inputLine = <ANNOINFILE>) {
        chomp $inputLine;
#        my ($chrom, $start, $stop, $ref, $var, $vcfChrom, $vcrPos, $vcfJunk, $vcrRef, $vcfVar, @junk);
        my ($chrom, $start, $stop, $ref, $var, $vcfJunk1, $vcfJunk2, $vcfJunk3, $vcfChrom, $vcrPos, $vcfJunk4, $vcrRef, $vcfVar, @junk);
        if ($RESUMEFOLDER && -e "$ANNOINFILE") {
            ($chrom, $start, $stop, $ref, $var, $vcfChrom, $vcrPos, $vcrRef, $vcfVar) = split /\t+/, $inputLine;
        } else {
#            ($chrom, $start, $stop, $ref, $var, $vcfChrom, $vcrPos, $vcfJunk, $vcrRef, $vcfVar, @junk) = split /\t+/, $inputLine;
            ($chrom, $start, $stop, $ref, $var, $vcfJunk1, $vcfJunk2, $vcfJunk3, $vcfChrom, $vcrPos, $vcfJunk4, $vcrRef, $vcfVar, @junk) = split /\t+/, $inputLine;
        }

       next if ($vcfVar eq ".");

        my $variantKey = $chrom . ":" . $start . "-" . $stop . "_" . $ref . "/" . $var; # " just to avoid a bug in eclipse display...

        # MULTIPLE CASES if start=stop then
        if ($start == $stop) {

            # INSERTION
            if ($ref eq '-') {
                $input{$variantKey}{VC} = 'insertion';

            # DELETION
            } elsif ($var eq '-') {
                $input{$variantKey}{VC} = 'deletion';

            # SUBSTITUTION
            } elsif (length ($var) > 1) {
                $input{$variantKey}{VC} = 'substitution';

            # SNP
            } else {
                $input{$variantKey}{VC} = 'SNP';
            }

        # DELETION if var='-'
        } elsif ($var eq '-') {
            $input{$variantKey}{VC} = 'deletion';

        # if SUBSTITUTION
        } else {
            $input{$variantKey}{VC} = 'substitution';
        }

        push @{$vcfVariants{$vcfChrom . ":" . $vcrPos . "_" . $vcrRef . "/" . $vcfVar}}, $variantKey;
        push @annoinfile, $chrom . "\t" . $start . "\t" . $stop . "\t" . $ref . "\t" . $var . "\t" . $vcfChrom . "\t" . $vcrPos . "\t" . $vcrRef . "\t" . $vcfVar . "\n";
    }
    close(ANNOINFILE);

    if (!($RESUMEFOLDER && -e "$ANNOINFILE")) {
        open(ANNOINFILE, ">$ANNOINFILE") || die "cannot open file $ANNOINFILE:\n$!\n";
        print ANNOINFILE @annoinfile;
        close(ANNOINFILE);
    }

    return (\%input, \%vcfVariants);
}

sub launchAnnovar {
    my %arg = @_;
    my $rO_vcf  = $arg{-vcf};
    my $rH_data = $arg{-data};

    $OUTFILE_ANNO .= "." . $BUILDVER;

    $USEDBUILD = $BUILDVER;
    ($USEDBUILD =~ /b37|v37/) and $USEDBUILD = "hg19";

    if (doIt(-elem => "gene")) {
        my $dbtype = ($MITO) ? "-dbtype MT_GRCh37_ensGene" : "";
        if ($RESUMEFOLDER && -e "$OUTFILE_ANNO.variant_function" && -e "$OUTFILE_ANNO.exonic_variant_function") {
            # resuming ANNOVAR -geneanno
            $VERBOSE1 and printerr("\n**********\n* Resuming gene annotation using already existing files $OUTFILE_ANNO.variant_function & $OUTFILE_ANNO.exonic_variant_function\n**********\n\n");
        } else {
            # run ANNOVAR -geneanno, first pass
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -geneanno $dbtype -buildver $USEDBUILD -neargene 10000 -splicing_threshold 6 -outfile $OUTFILE_ANNO $ANNOINFILE $GENANNOPATH -separate -exonicsplicing -transcript_function $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting genome annotation with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Gene annotation done...\n**********\n\n");
        }

        # open ANNOVAR output files for processing
        open(ANNOVAR_MAIN_GENERAL, "<$OUTFILE_ANNO.variant_function") || die "cannot open file $OUTFILE_ANNO.variant_function:\n$!\n";
        open(ANNOVAR_MAIN_DETAILED, "<$OUTFILE_ANNO.exonic_variant_function") || die "cannot open file $OUTFILE_ANNO.exonic_variant_function:\n$!\n";

        $rO_vcf->add_header_line({key=>'INFO', ID=>'VC',   Number=>'A', Type=>'String', Description=>'Variant class'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'VFT',  Number=>'A', Type=>'String', Description=>'Variant function type'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'VT',   Number=>'A', Type=>'String', Description=>'Coding variant type'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'Gene', Number=>1,   Type=>'String', Description=>'Gene symbol'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'DA',   Number=>'A', Type=>'String', Description=>'Detailed annotation of the variant'});

        # annotate noncoding variants - get gene symbol (col 2) and general functional effect (col 1)
        $VERBOSE1 and printerr("\n**********\n* Parsing $OUTFILE_ANNO.variant_function to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_MAIN_GENERAL>) {
            if ($USEDBUILD =~ /mm\d+/) {
                chomp $line;  # $1     $2      $3   $4         $5     $6       $7          $8
                if ($line =~ /^(.+?)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+).+/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    my $varFuncType = $1;
                    my $gene = $2;
                    $gene =~ s/;/:/g;   # ABC;NM_blabla becomes ABC:NM_blabla
                    push @{$rH_data->{$key}->{'VFT'}}, $varFuncType || "";
                    if ($varFuncType eq "intergenic") {
                        $gene =~ s/dist=/dist:/g;
                        $gene =~ s/,/-/g;
                        $gene =~ s/(\w+)\:\w+(\(dist:\d+\))/$1$2/g;
                        push @{$rH_data->{$key}->{'DA'}}, "$varFuncType:$gene";
                    } else {
                        $gene =~ s/,/|/g;
                        push @{$rH_data->{$key}->{'Gene'}}, $gene;
                        my @genes = split(/\|/, $gene);
                        $gene = join "|$varFuncType:", sort @genes;
                        if ($varFuncType !~ /^exonic/) {
                            push @{$rH_data->{$key}->{'DA'}}, "$varFuncType:$gene";
                        }
                    }
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            } else {
                chomp $line;  # $1     $2     $3   $4          $5     $6       $7          $8
                if ($line =~ /^(.+?)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    my $varFuncType = $1;
                    my $gene = $2;
                    $gene =~ s/;/:/g;	# ABC;NM_blabla becomes ABC:NM_blabla
                    push @{$rH_data->{$key}->{'VFT'}}, $varFuncType || "";
                    if ($varFuncType eq "intergenic") {
                        $gene =~ s/dist=/dist:/g;
                        $gene =~ s/,/-/g;
                        $gene =~ s/(\w+)\:\w+(\(dist:\d+\))/$1$2/g;
                        push @{$rH_data->{$key}->{'DA'}}, "$varFuncType:$gene";
                    } else {
                        $gene =~ s/,/|/g;
                        push @{$rH_data->{$key}->{'Gene'}}, $gene;
                        my @genes = split /\|/, $gene;  # / just to avoid a bug in eclipse display...
                        $gene = join "|$varFuncType:", sort @genes;
                        if ($varFuncType !~ /^exonic/) {
                            push @{$rH_data->{$key}->{'DA'}}, "$varFuncType:$gene";
                        }
                    }
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
        }
        close(ANNOVAR_MAIN_GENERAL);

        # further annotate coding variants - using line number (col 1) to compare to key of global input hash, annotate synonymous/nonsynonymous/frameshift/nonframeshift/insertion/deletion (col 2) and isoform/protein/cdna/AA info (col 3)
        $VERBOSE1 and printerr("\n**********\n* Parsing $OUTFILE_ANNO.exonic_variant_function to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_MAIN_DETAILED>) {
            if ($USEDBUILD =~ /mm\d+/) {
                chomp $line;    #        $1     $2     $3    $4        $5     $6       $7          $8
                if ($line =~ /^line\d+\t(.+?)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+).+/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    my $variantType = $1;
                    my $annotation = $2;
                    $annotation =~ s/,\s*$//g;
                    $annotation =~ s/,/||/g;
                    $variantType =~ s/ /_/;
                    my $annotString = "";
                    foreach my $varFuncType (@{$rH_data->{$key}->{'VFT'}}) {
                        if ($varFuncType eq "exonic" || $varFuncType eq "exonic_splicing") {
                            my @annotations = split(/\|\|/, $annotation);
                            $annotation = join "|$varFuncType:$variantType:", sort @annotations;
                            $annotString = "$varFuncType:$variantType:$annotation";
                        }
                    }
                    if ($annotString eq "") {
                        my @annotations = split(/\|\|/, $annotation);
                        $annotation = join "|exonic:$variantType:", sort @annotations;
                    }
                    push @{$rH_data->{$key}->{'DA'}}, $annotString;
                    $rH_data->{$key}->{'VT'} = $variantType;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            } else {
                chomp $line;    #        $1     $2      $3    $4        $5     $6       $7         $8
                if ($line =~ /^line\d+\t(.+?)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    my $variantType = $1;
                    my $annotation = $2;
                    $annotation =~ s/,\s*$//g;
                    $annotation =~ s/,/||/g;
                    $variantType =~ s/ /_/;
                    my $annotString = "";
                    foreach my $varFuncType (@{$rH_data->{$key}->{'VFT'}}) {
                        if ($varFuncType eq "exonic" || $varFuncType eq "exonic_splicing") {
                            my @annotations = split /\|\|/, $annotation;    # / just to avoid a bug in eclipse display...
                            $annotation = join "|$varFuncType:$variantType:", sort @annotations;
                            $annotString = "$varFuncType:$variantType:$annotation";
                        }
                    }
                    if ($annotString eq "") {
                        my @annotations = split /\|\|/, $annotation;    # / just to avoid a bug in eclipse display...
                        $annotation = join "|exonic:$variantType:", sort @annotations;
                    }
                    push @{$rH_data->{$key}->{'DA'}}, $annotString;
                    $rH_data->{$key}->{'VT'} = $variantType;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
        }
        close(ANNOVAR_MAIN_DETAILED);
    }

    # if not excluded, run ANNOVAR -filter to output dbSNP132 variants
    if (doIt(-elem => "snp132")) {
        my $annovarOutput = "$OUTFILE_DBSNP132.$USEDBUILD\_snp132_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_132', Number=>1, Type=>'String', Description=>"All SNPs(132) - all SNPs from dbSNP mapping to reference assembly. (version snp132)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype snp132
            $VERBOSE1 and printerr("\n**********\n* Resuming filtering against dbSNP132 using already existing file $annovarOutput\n**********\n\n");
        } else {
            # First run annovar
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype snp132 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP132 $ANNOINFILE $DBSNPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against dbSNP132 with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against dbSNP132 done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_DBSNP132, "<$annovarOutput") || die "cannot open file $annovarOutput:\n$!\n";

        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_DBSNP132>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;  #   $1       $2      $3   $4        $5     $6       $7         $8
            if ($line =~ /(snp132)\t(rs\d+)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{dbSNP_132} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSNP132);
    }


    # if not excluded, run ANNOVAR -filter to output dbSNP138 variants
    if (doIt(-elem => "snp138")) {
        my $annovarOutput = "$OUTFILE_DBSNP138.$USEDBUILD\_snp138_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_138', Number=>1, Type=>'String', Description=>"All SNPs(138) - all SNPs from dbSNP mapping to reference assembly. (version snp138)"});
        if ($RESUMEFOLDER && -e "$annovarOutput") {
            # resuming ANNOVAR -filter -dbtype snp138
            $VERBOSE1 and printerr("\n**********\n* Resuming filtering against dbSNP138 using already existing file $annovarOutput\n**********\n\n");
        } else {
            # First run annovar 
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype snp138 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP138 $ANNOINFILE $DBSNPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against dbSNP138 with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against dbSNP138 done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_DBSNP138, "<$annovarOutput") || die "cannot open file $annovarOutput:\n$!\n";

        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_DBSNP138>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #     $1    2     $3    $4         $5     $6       $7         $8
            if ($line =~ /(snp138)\t(rs\d+)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{dbSNP_138} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSNP138);
    }

    # if not excluded, run ANNOVAR -filter to output dbSNP138Common variants
    if (doIt(-elem => "snp138common")) {
        my $dbFile = "hg19_snp138Common.txt.generic";
        my $annovarOutput = "$OUTFILE_DBSNP138_COMMON.$USEDBUILD\_generic_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_138Common', Number=>0, Type=>'Flag', Description=>'Common SNPs(138) - SNPs with more than 1% minor allele frequency (MAF), mapping only once to reference assembly. (version snp138)'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_138CommonAlCount', Number=>1, Type=>'Float', Description=>"All Allele Count for Common SNPs(138) (version snp138Common)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_138CommonAlFreq', Number=>1, Type=>'Float', Description=>"Mutant Allele Frequency for Common SNPs(138) (version snp138Common)"});
        if ($RESUMEFOLDER && -e "$annovarOutput") {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming filtering against dbSNP_138Common using already existing file $annovarOutput\n**********\n\n");
        } else {
            # First run annovar
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbFile -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP138_COMMON $ANNOINFILE $DBSNPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against dbSNP_138Common with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against dbSNP_138Common done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_COMMONDBSNP, "<$annovarOutput") || die "cannot open file $annovarOutput:\n$!\n";

        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_COMMONDBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1      $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t(rs\d+.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($rs, $aCount, $aFreq) = split(/,/, $2);
                if (defined $rH_data->{$key}->{dbSNP_138} && $rH_data->{$key}->{dbSNP_138} =~ /^rs.+/) {
                    (defined $rs) and $rH_data->{$key}->{dbSNP_138Common} = 1;
                    (defined $aCount) and $rH_data->{$key}->{dbSNP_138CommonAlCount} = $aCount;
                    (defined $aFreq) and $rH_data->{$key}->{dbSNP_138CommonAlFreq} = $aFreq;
                }
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_COMMONDBSNP);
    }

    # if not excluded, run ANNOVAR -filter to output dbSNP138Common variants
    if (doIt(-elem => "snp138flagged")) {
        my $dbFile = "hg19_snp138Flagged.txt.generic";
        my $annovarOutput = "$OUTFILE_DBSNP138_FLAGGED.$USEDBUILD\_generic_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_138Flagged', Number=>0, Type=>'Flag', Description=>"Flagged SNPs(138) - SNPs less than 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as 'clinically associated' (version snp138)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_138FlaggedAlCount', Number=>1, Type=>'Float', Description=>"All Allele Count for Flagged SNPs(138) (version snp138Flagged)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_138FlaggedAlFreq', Number=>1, Type=>'Float', Description=>"Mutant Allele Frequency for Flagged SNPs(138) (version snp138Flagged)"});
        if ($RESUMEFOLDER && -e "$annovarOutput") {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming filtering against dbSNP138_Flagged using already existing file $annovarOutput\n**********\n\n");
        } else {
            # First run annovar
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbFile -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP138_FLAGGED $ANNOINFILE $DBSNPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against dbSNP_138Flagged with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against dbSNP138_Flagged done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_CLINICDBSNP, "<$annovarOutput") || die "cannot open file $annovarOutput:\n$!\n";

        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_CLINICDBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1         $2      $3     $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t(rs\d+.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($rs, $aCount, $aFreq) = split(/,/, $2);
                if (defined $rH_data->{$key}->{dbSNP_138} && $rH_data->{$key}->{dbSNP_138} =~ /^rs.+/) {
                    (defined $rs) and $rH_data->{$key}->{dbSNP_138Flagged} = 1;
                    (defined $aCount) and $rH_data->{$key}->{dbSNP_138FlaggedAlCount} = $aCount;
                    (defined $aFreq) and $rH_data->{$key}->{dbSNP_138FlaggedAlFreq} = $aFreq;
                }
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_CLINICDBSNP);
    }


    # if not excluded, run ANNOVAR -filter to output dbSNP142 variants
    if (doIt(-elem => "snp142")) {
        my $annovarOutput = "$OUTFILE_DBSNP142.$USEDBUILD\_snp142_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_142', Number=>1, Type=>'String', Description=>"All SNPs(142) - all SNPs from dbSNP mapping to reference assembly. (version snp142)"});
        if ($RESUMEFOLDER && -e "$annovarOutput") {
            # resuming ANNOVAR -filter -dbtype snp142
            $VERBOSE1 and printerr("\n**********\n* Resuming filtering against dbSNP142 using already existing file $annovarOutput\n**********\n\n");
        } else {
            # First run annovar 
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype snp142 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP142 $ANNOINFILE $DBSNPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against dbSNP142 with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against dbSNP142 done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_DBSNP142, "<$annovarOutput") || die "cannot open file $annovarOutput:\n$!\n";

        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_DBSNP142>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #     $1    2     $3    $4         $5     $6       $7         $8
            if ($line =~ /(snp142)\t(rs\d+)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{dbSNP_142} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSNP142);
    }

    # if not excluded, run ANNOVAR -filter to output dbSNP142Common variants
    if (doIt(-elem => "snp142common")) {
        my $dbFile = "hg19_snp142Common.txt.generic";
        my $annovarOutput = "$OUTFILE_DBSNP142_COMMON.$USEDBUILD\_generic_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_142Common', Number=>0, Type=>'Flag', Description=>'Common SNPs(142) - SNPs with more than 1% minor allele frequency (MAF), mapping only once to reference assembly. (version snp142)'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_142CommonAlCount', Number=>1, Type=>'Float', Description=>"All Allele Count for Common SNPs(142) (version snp142Common)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_142CommonAlFreq', Number=>1, Type=>'Float', Description=>"Mutant Allele Frequency for Common SNPs(142) (version snp142Common)"});
        if ($RESUMEFOLDER && -e "$annovarOutput") {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming filtering against dbSNP_142Common using already existing file $annovarOutput\n**********\n\n");
        } else {
            # First run annovar
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbFile -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP142_COMMON $ANNOINFILE $DBSNPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against dbSNP_142Common with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against dbSNP_142Common done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_COMMONDBSNP, "<$annovarOutput") || die "cannot open file $annovarOutput:\n$!\n";

        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_COMMONDBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1      $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t(rs\d+.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($rs, $aCount, $aFreq) = split(/,/, $2);
                if (defined $rH_data->{$key}->{dbSNP_142} && $rH_data->{$key}->{dbSNP_142} =~ /^rs.+/) {
                    (defined $rs) and $rH_data->{$key}->{dbSNP_142Common} = 1;
                    (defined $aCount) and $rH_data->{$key}->{dbSNP_142CommonAlCount} = $aCount;
                    (defined $aFreq) and $rH_data->{$key}->{dbSNP_142CommonAlFreq} = $aFreq;
                }
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_COMMONDBSNP);
    }

    # if not excluded, run ANNOVAR -filter to output dbSNP142Common variants
    if (doIt(-elem => "snp142flagged")) {
        my $dbFile = "hg19_snp142Flagged.txt.generic";
        my $annovarOutput = "$OUTFILE_DBSNP142_FLAGGED.$USEDBUILD\_generic_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_142Flagged', Number=>0, Type=>'Flag', Description=>"Flagged SNPs(142) - SNPs less than 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as 'clinically associated' (version snp142)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_142FlaggedAlCount', Number=>1, Type=>'Float', Description=>"All Allele Count for Flagged SNPs(142) (version snp142Flagged)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbSNP_142FlaggedAlFreq', Number=>1, Type=>'Float', Description=>"Mutant Allele Frequency for Flagged SNPs(142) (version snp142Flagged)"});
        if ($RESUMEFOLDER && -e "$annovarOutput") {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming filtering against dbSNP142_Flagged using already existing file $annovarOutput\n**********\n\n");
        } else {
            # First run annovar
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbFile -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP142_FLAGGED $ANNOINFILE $DBSNPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against dbSNP_142Flagged with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against dbSNP142_Flagged done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_CLINICDBSNP, "<$annovarOutput") || die "cannot open file $annovarOutput:\n$!\n";

        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_CLINICDBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1         $2      $3     $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t(rs\d+.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($rs, $aCount, $aFreq) = split(/,/, $2);
                if (defined $rH_data->{$key}->{dbSNP_142} && $rH_data->{$key}->{dbSNP_142} =~ /^rs.+/) {
                    (defined $rs) and $rH_data->{$key}->{dbSNP_142Flagged} = 1;
                    (defined $aCount) and $rH_data->{$key}->{dbSNP_142FlaggedAlCount} = $aCount;
                    (defined $aFreq) and $rH_data->{$key}->{dbSNP_142FlaggedAlFreq} = $aFreq;
                }
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_CLINICDBSNP);
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg")) {
        my $dbtype = "1000g2015aug_all";
        my $annovarOutput = "$OUTFILE_1000G.$USEDBUILD\_ALL.sites.2015_08_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KG', Number=>'A', Type=>'Float', Description=>"Frequency in the 1000 Genome project (201508 collection v5b (based on 201305 alignment) with including chrX and chrY data)"});
                                                                                                                                #         Autosomes from 2013_05_02 version & chromosome X from 2012.04 collection v4 (based on 2011.03 alignment)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype 1000g2015aug_all
            $VERBOSE1 and printerr("\n**********\n* Resuming 1000 Genomes frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_1000G $ANNOINFILE $FREQUENCIESPATH/1000g/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against 1000 Genomes db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against 1000 Genomes db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_1000g, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_1000g>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1          $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t([01]\.*\d*)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{KG} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_1000g);
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg_afr")) {
        my $dbtype = "1000g2015aug_afr";
        my $annovarOutput = "$OUTFILE_AFR.$USEDBUILD\_AFR.sites.2015_08_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KG_AFR', Number=>'A', Type=>'Float', Description=>"Allele frequency of African population of the 1000 Genome project (201508 collection v5b (based on 201305 alignment) with including chrX and chrY data)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype 1000g2015aug_afr
            $VERBOSE1 and printerr("\n**********\n* Resuming 1000 Genomes African frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AFR $ANNOINFILE $FREQUENCIESPATH/1000g/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against 1000 Genomes African db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against 1000 Genomes African frequencies db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_AFR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_AFR>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1          $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t([01]\.*\d*)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{KG_AFR} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_AFR);
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg_amr")) {
        my $dbtype = "1000g2015aug_amr";
        my $annovarOutput = "$OUTFILE_AMR.$USEDBUILD\_AMR.sites.2015_08_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KG_AMR', Number=>'A', Type=>'Float', Description=>"Allele frequency of Ad Mixed American population of the 1000 Genome project (201508 collection v5b (based on 201305 alignment) with including chrX and chrY data)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype 1000g2015aug_amr
            $VERBOSE1 and printerr("\n**********\n* Resuming 1000 Genomes Ad Mixed American frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AMR $ANNOINFILE $FREQUENCIESPATH/1000g/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against 1000 Genomes Ad Mixed American db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against 1000 Genomes Ad Mixed American frequencies db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_AMR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_AMR>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1          $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t([01]\.*\d*)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{KG_AMR} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_AMR);
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg_eur")) {
        my $dbtype = "1000g2015aug_eur";
        my $annovarOutput = "$OUTFILE_EUR.$USEDBUILD\_EUR.sites.2015_08_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KG_EUR', Number=>'A', Type=>'Float', Description=>"Allele frequency of European population of the 1000 Genome project (201508 collection v5b (based on 201305 alignment) with including chrX and chrY data)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype 1000g2015aug_eur
            $VERBOSE1 and printerr("\n**********\n* Resuming 1000 Genomes European frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EUR $ANNOINFILE $FREQUENCIESPATH/1000g/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against 1000 Genomes European db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against 1000 Genomes European frequencies db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_EUR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_EUR>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1          $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t([01]\.*\d*)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{KG_EUR} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_EUR);
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg_eas")) {
        my $dbtype = "1000g2015aug_eas";
        my $annovarOutput = "$OUTFILE_EAS.$USEDBUILD\_EAS.sites.2015_08_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KG_EAS', Number=>'A', Type=>'Float', Description=>"Allele frequency of East Asian population of the 1000 Genome project (201508 collection v5b (based on 201305 alignment) with including chrX and chrY data)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype 1000g2015aug_eas
            $VERBOSE1 and printerr("\n**********\n* Resuming 1000 Genomes East Asian frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EAS $ANNOINFILE $FREQUENCIESPATH/1000g/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against 1000 Genomes East Asian db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against 1000 Genomes East Asian frequencies db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_EAS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_EAS>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1          $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t([01]\.*\d*)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{KG_EAS} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_EAS);
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg_sas")) {
        my $dbtype = "1000g2015aug_sas";
        my $annovarOutput = "$OUTFILE_SAS.$USEDBUILD\_SAS.sites.2015_08_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KG_SAS', Number=>'A', Type=>'Float', Description=>"Allele frequency of South Asian population of the 1000 Genome project (201508 collection v5b (based on 201305 alignment) with including chrX and chrY data)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype 1000g2015aug_sas
            $VERBOSE1 and printerr("\n**********\n* Resuming 1000 Genomes South Asian frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_SAS $ANNOINFILE $FREQUENCIESPATH/1000g/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against 1000 Genomes South Asian db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against 1000 Genomes South Asian frequencies db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_SAS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_SAS>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1          $2        $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t([01]\.*\d*)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{KG_SAS} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_SAS);
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "exac")) {
        my $dbtype = "exac03";
        my $annovarOutput = "$OUTFILE_EXAC.$USEDBUILD\_exac03_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_Freq', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data (version 0.3 - 2015.07.29)"});
        (doIt(-elem => "exac_afr")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_AFR', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data in African population (version 0.3 - 2015.07.29.)"});
        (doIt(-elem => "exac_amr")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_AMR', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data in Admixed American population (version 0.3 - 2015.07.29.)"});
        (doIt(-elem => "exac_eas")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_EAS', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data in East Asian population (version 0.3 - 2015.07.29.)"});
        (doIt(-elem => "exac_fin")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_FIN', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data in Finnish population (version 0.3 - 2015.07.29.)"});
        (doIt(-elem => "exac_nfe")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_NFE', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data in Non-finnish European population (version 0.3 - 2015.07.29.)"});
        (doIt(-elem => "exac_oth")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_OTH', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data in Other population (version 0.3 - 2015.07.29.)"});
        (doIt(-elem => "exac_sas")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_SAS', Number=>'A', Type=>'Float', Description=>"ExAC 65000 exome allele frequency data in South Asian population (version 0.3 - 2015.07.29.)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype $dbtype
            $VERBOSE1 and printerr("\n**********\n* Resuming ExAC frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EXAC $ANNOINFILE $FREQUENCIESPATH/ExAC/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against ExAC db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against ExAC db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_EXAC, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_EXAC>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1       $2    $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($ExAC_Freq, $ExAC_AFR, $ExAC_AMR, $ExAC_EAS, $ExAC_FIN, $ExAC_NFE, $ExAC_OTH, $ExAC_SAS) = split /,/, $2;
                $rH_data->{$key}->{ExAC_Freq} = ($ExAC_Freq eq ".") ? undef : $ExAC_Freq;
                (doIt(-elem => "exac_afr")) and $rH_data->{$key}->{ExAC_AFR}  = ($ExAC_AFR eq ".")  ? undef : $ExAC_AFR;
                (doIt(-elem => "exac_amr")) and $rH_data->{$key}->{ExAC_AMR}  = ($ExAC_AMR eq ".")  ? undef : $ExAC_AMR;
                (doIt(-elem => "exac_eas")) and $rH_data->{$key}->{ExAC_EAS}  = ($ExAC_EAS eq ".")  ? undef : $ExAC_EAS;
                (doIt(-elem => "exac_fin")) and $rH_data->{$key}->{ExAC_FIN}  = ($ExAC_FIN eq ".")  ? undef : $ExAC_FIN;
                (doIt(-elem => "exac_nfe")) and $rH_data->{$key}->{ExAC_NFE}  = ($ExAC_NFE eq ".")  ? undef : $ExAC_NFE;
                (doIt(-elem => "exac_oth")) and $rH_data->{$key}->{ExAC_OTH}  = ($ExAC_OTH eq ".")  ? undef : $ExAC_OTH;
                (doIt(-elem => "exac_sas")) and $rH_data->{$key}->{ExAC_SAS}  = ($ExAC_SAS eq ".")  ? undef : $ExAC_SAS;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_EXAC);
    }
    # if not excluded, run ANNOVAR -filter to output gme variants
    if (doIt(-elem => "gme")) {
        my $dbtype = "gme";
        my $annovarOutput = "$OUTFILE_GME.$USEDBUILD\_gme_dropped";

        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_AF',    Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency (version Annovar - 20161024)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_NWA',   Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency in NorthWest Africa (version Annovar - 20161024)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_NEA',   Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency in NorthEast Africa (version Annovar - 20161024)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_AP',    Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency in Arabian peninsula (version Annovar - 20161024)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_Israel',Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency in Israel (version Annovar - 20161024)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_SD',    Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency in Syrian desert (version Annovar - 20161024)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_TP',    Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency in Turkish peninsula (version Annovar - 20161024)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GME_CA',    Number=>'A', Type=>'Float', Description=>"Great Middle East allele frequency in Central Asia (version Annovar - 20161024)"});
        if ($RESUMEFOLDER && -e $annovarOutput) {
          # resuming ANNOVAR -filter -dbtype $dbtype
            $VERBOSE1 and printerr("\n**********\n* Resuming GME frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GME $ANNOINFILE $GMEPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against GME db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against GME db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
       open(ANNOVAR_GME, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GME>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1       $2    $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($GME_AF, $GME_NWA, $GME_NEA, $GME_AP, $GME_Israel, $GME_SD, $GME_TP, $GME_CA) = split /,/, $2;
                $rH_data->{$key}->{GME_AF}   = ($GME_AF eq ".") ? undef : $GME_AF;
                $rH_data->{$key}->{GME_NWA}  = ($GME_NWA eq ".")  ? undef : $GME_NWA;
                $rH_data->{$key}->{GME_NEA}  = ($GME_NEA eq ".")  ? undef : $GME_NEA;
                $rH_data->{$key}->{GME_AP}  = ($GME_AP eq ".")  ? undef : $GME_AP;
                $rH_data->{$key}->{GME_Israel}  = ($GME_Israel eq ".")  ? undef : $GME_Israel;
                $rH_data->{$key}->{GME_SD}  = ($GME_SD eq ".")  ? undef : $GME_SD;
                $rH_data->{$key}->{GME_TP}  = ($GME_TP eq ".")  ? undef : $GME_TP;
                $rH_data->{$key}->{GME_CA}  = ($GME_CA eq ".")  ? undef : $GME_CA;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GME);
    }                             

    # if not excluded, run ANNOVAR -filter to output gnomad exome variants
    if (doIt(-elem => "gnomad_exome")) {
        my $dbtype ="gnomad_exome";
        my $annovarOutput = "$OUTFILE_GNOMAD_EXOME.$USEDBUILD\_gnomad_exome_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_ALL',  Number=>'A', Type=>'Float', Description=>"gnomAD ALL exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_AFR',  Number=>'A', Type=>'Float', Description=>"gnomAD AFR exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_AMR',  Number=>'A', Type=>'Float', Description=>"gnomAD AMR exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_ASJ',  Number=>'A', Type=>'Float', Description=>"gnomAD ASJ exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_EAS',  Number=>'A', Type=>'Float', Description=>"gnomAD EAS exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_FIN',  Number=>'A', Type=>'Float', Description=>"gnomAD FIN exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_NFE',  Number=>'A', Type=>'Float', Description=>"gnomAD NFE exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_OTH',  Number=>'A', Type=>'Float', Description=>"gnomAD OTH exome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_exome_SAS',  Number=>'A', Type=>'Float', Description=>"gnomAD SAS exome allele frequency (version Annovar - 20170311)"});

        if ($RESUMEFOLDER && -e $annovarOutput) {
           # resuming ANNOVAR -filter -dbtype $dbtype
            $VERBOSE1 and printerr("\n**********\n* Resuming gnomAD exome frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GNOMAD_EXOME $ANNOINFILE $GNOMADPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against gnomAD exome db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against gnomAD exome db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
       open(ANNOVAR_GNOMAD_EXOME, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GNOMAD_EXOME>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1       $2    $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($gnomAD_exome_ALL, $gnomAD_exome_AFR, $gnomAD_exome_AMR, $gnomAD_exome_ASJ, $gnomAD_exome_EAS, $gnomAD_exome_FIN, $gnomAD_exome_NFE, $gnomAD_exome_OTH, $gnomAD_exome_SAS) = split /,/, $2;
                $rH_data->{$key}->{gnomAD_exome_ALL}  = ($gnomAD_exome_ALL eq ".") ? undef : $gnomAD_exome_ALL;
                $rH_data->{$key}->{gnomAD_exome_AFR}  = ($gnomAD_exome_AFR eq ".")  ? undef : $gnomAD_exome_AFR;
                $rH_data->{$key}->{gnomAD_exome_AMR}  = ($gnomAD_exome_AMR eq ".")  ? undef : $gnomAD_exome_AMR;
                $rH_data->{$key}->{gnomAD_exome_ASJ}  = ($gnomAD_exome_ASJ eq ".")  ? undef : $gnomAD_exome_ASJ;
                $rH_data->{$key}->{gnomAD_exome_EAS}  = ($gnomAD_exome_EAS eq ".")  ? undef : $gnomAD_exome_EAS;
                $rH_data->{$key}->{gnomAD_exome_FIN}  = ($gnomAD_exome_FIN eq ".")  ? undef : $gnomAD_exome_FIN;
                $rH_data->{$key}->{gnomAD_exome_NFE}  = ($gnomAD_exome_NFE eq ".")  ? undef : $gnomAD_exome_NFE;
                $rH_data->{$key}->{gnomAD_exome_OTH}  = ($gnomAD_exome_OTH eq ".")  ? undef : $gnomAD_exome_OTH;
                $rH_data->{$key}->{gnomAD_exome_SAS}  = ($gnomAD_exome_SAS eq ".")  ? undef : $gnomAD_exome_SAS;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GNOMAD_EXOME);
    }
                               
    # if not excluded, run ANNOVAR -filter to output gnomad wgs variants
    if (doIt(-elem => "gnomad_genome")) {
        my $dbtype ="gnomad_genome";
        my $annovarOutput = "$OUTFILE_GNOMAD_WGS.$USEDBUILD\_$dbtype\_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_ALL',  Number=>'A', Type=>'Float', Description=>"gnomAD ALL genome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_AFR',  Number=>'A', Type=>'Float', Description=>"gnomAD AFR genome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_AMR',  Number=>'A', Type=>'Float', Description=>"gnomAD AMR genome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_ASJ',  Number=>'A', Type=>'Float', Description=>"gnomAD ASJ genome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_EAS',  Number=>'A', Type=>'Float', Description=>"gnomAD EAS genome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_FIN',  Number=>'A', Type=>'Float', Description=>"gnomAD FIN genome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_NFE',  Number=>'A', Type=>'Float', Description=>"gnomAD NFE genome allele frequency (version Annovar - 20170311)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'gnomAD_genome_OTH',  Number=>'A', Type=>'Float', Description=>"gnomAD OTH genome allele frequency (version Annovar - 20170311)"});

        if ($RESUMEFOLDER && -e $annovarOutput) {
           # resuming ANNOVAR -filter -dbtype $dbtype
            $VERBOSE1 and printerr("\n**********\n* Resuming gnomAD genome frequencies annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GNOMAD_WGS $ANNOINFILE $GNOMADPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against gnomAD genome db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against gnomAD genome db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
       open(ANNOVAR_GNOMAD_WGS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GNOMAD_WGS>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1       $2    $3    $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($gnomAD_genome_ALL, $gnomAD_genome_AFR, $gnomAD_genome_AMR, $gnomAD_genome_ASJ, $gnomAD_genome_EAS, $gnomAD_genome_FIN, $gnomAD_genome_NFE, $gnomAD_genome_OTH) = split /,/, $2;
                $rH_data->{$key}->{gnomAD_genome_ALL}  = ($gnomAD_genome_ALL eq ".") ? undef : $gnomAD_genome_ALL;
                $rH_data->{$key}->{gnomAD_genome_AFR}  = ($gnomAD_genome_AFR eq ".")  ? undef : $gnomAD_genome_AFR;
                $rH_data->{$key}->{gnomAD_genome_AMR}  = ($gnomAD_genome_AMR eq ".")  ? undef : $gnomAD_genome_AMR;
                $rH_data->{$key}->{gnomAD_genome_ASJ}  = ($gnomAD_genome_ASJ eq ".")  ? undef : $gnomAD_genome_ASJ;
                $rH_data->{$key}->{gnomAD_genome_EAS}  = ($gnomAD_genome_EAS eq ".")  ? undef : $gnomAD_genome_EAS;
                $rH_data->{$key}->{gnomAD_genome_FIN}  = ($gnomAD_genome_FIN eq ".")  ? undef : $gnomAD_genome_FIN;
                $rH_data->{$key}->{gnomAD_genome_NFE}  = ($gnomAD_genome_NFE eq ".")  ? undef : $gnomAD_genome_NFE;
                $rH_data->{$key}->{gnomAD_genome_OTH}  = ($gnomAD_genome_OTH eq ".")  ? undef : $gnomAD_genome_OTH;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GNOMAD_WGS);
    }
                               

    # if not excluded, run ANNOVAR -filter to output Complete Genomics variants
    if (doIt(-elem => "cg")) {
        my $annovarOutput = "$OUTFILE_CG46.$USEDBUILD\_cg46\_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'CG46', Number=>'A', Type=>'Float', Description=>'allele frequency in 46 unrelated human subjects sequenced by Complete Genomics (db updated on 20140106)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype cg46
            $VERBOSE1 and printerr("\n**********\n* Resuming Complete Genomics annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype cg46 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_CG46 $ANNOINFILE $CGPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting filtering against Complete Genomics db with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Filtering against Complete Genomics db done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_CG, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_CG>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1         $2     $3    $4        $5     $6       $7         $8
            if ($line =~ /(cg46)\t([01]\.\d*)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{CG46} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_CG);
    }

    # if not excluded, run ANNOVAR -filter to output COSMIC annotations
    if (doIt(-elem => "cosmic")) {
        my $annovarOutput = "$OUTFILE_COSMIC.$USEDBUILD\_cosmic70_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'COSMIC', Number=>1, Type=>'String', Description=>'Sanger COSMIC (Catalogue of Somatic Mutations in Cancer) version 70 - somatic disease mutations description (ANNOVAR - cosmic70)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype cosmic70
            $VERBOSE1 and printerr("\n**********\n* Resuming COSMIC annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype cosmic70 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_COSMIC $ANNOINFILE $CANCERPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of COSMIC annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of COSMIC annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_COSMIC, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_COSMIC>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(cosmic70)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{COSMIC} = $2;
                $rH_data->{$key}->{COSMIC} =~ s/\=/\:/g;
                $rH_data->{$key}->{COSMIC} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{COSMIC} =~ s/ /_/g;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_COSMIC);
    }

    # if not excluded, run ANNOVAR -filter to output ICGC21 annotations
    if (doIt(-elem => "icgc21")) {
        my $annovarOutput = "$OUTFILE_ICGC21.$USEDBUILD\_icgc21_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ICGC21_ID', Number=>1, Type=>'String', Description=>'International Cancer Genome Consortium version 21 - ID (ANNOVAR - 20160622)'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ICGC21_Occurrence', Number=>1, Type=>'String', Description=>'International Cancer Genome Consortium version 21 - Occurrence (ANNOVAR - 20160622)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype icgc21    
            $VERBOSE1 and printerr("\n**********\n* Resuming ICGC21 annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype icgc21 -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ICGC21 $ANNOINFILE $CANCERPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ICGC21 annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of ICGC21 annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ICGC21, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ICGC21>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(icgc21)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($icgc_id, $icgc_occurrence) = split(/,/, $2);
                $rH_data->{$key}->{ICGC21_ID} = $icgc_id;
                $rH_data->{$key}->{ICGC21_Occurrence} = $icgc_occurrence;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ICGC21);
    }
                


    # if not excluded, run ANNOVAR -filter to output WGSDANN annotations
    if (doIt(-elem => "wgsdann")) {
        my $annovarOutput = "$OUTFILE_WGSDANN.$USEDBUILD\_dann_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'WGSDANN', Number=>1, Type=>'Float', Description=>'WGSDANN scores whole-genome variants by training a deep neural network (DNN). DNNs can capture non-linear relationships among features and are better suited than SVMs for problems with a large number of samples and features - version_201410'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype dann
            $VERBOSE1 and printerr("\n**********\n* Resuming WGSDANN annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype dann -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_WGSDANN $ANNOINFILE $WGSDANNPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of WGSDANN annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of WGSDANN annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_WGSDANN, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
	$VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_WGSDANN>) {         
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(dann)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{WGSDANN} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_WGSDANN);
    }

    # if not excluded, run ANNOVAR -filter to output GWAVA annotations
    if (doIt(-elem => "gwava")) {
        my $annovarOutput = "$OUTFILE_GWAVA.$USEDBUILD\_gwava_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GWAVA_region_score', Number=>1, Type=>'Float', Description=>'whole genome GWAVA_region_score - version_20150623'});
	$rO_vcf->add_header_line({key=>'INFO', ID=>'GWAVA_tss_score', Number=>1, Type=>'Float', Description=>'whole genome GWAVA_tss_score - version_20150623'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype gwava
            $VERBOSE1 and printerr("\n**********\n* Resuming GWAVA annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype gwava -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GWAVA $ANNOINFILE $GWAVAPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GWAVA annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GWAVA annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GWAVA, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GWAVA>) {
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(gwava)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
		my ($region, $tss, $unmatch) = split(/,/, $2);
                $rH_data->{$key}->{GWAVA_region_score} = $region;
		$rH_data->{$key}->{GWAVA_tss_score} = $tss;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GWAVA);
    }


    # if not excluded, run ANNOVAR -filter to output WGSCADD annotations
    if (doIt(-elem => "wgscadd")) {
        my $annovarOutput = "$OUTFILE_WGSCADD.$USEDBUILD\_cadd13_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'WGSCADD', Number=>1, Type=>'Float', Description=>'CADD 1.3 Raw Score (Combined Annotation Dependent Depletion) is a score that is based on SVM on multiple other scores. It assigns a score to each possible mutation in the human genome, therefore can evaluate non-coding variants as well as coding ones - version 1.3 - Annovar 20170123'});
	$rO_vcf->add_header_line({key=>'INFO', ID=>'WGSCADD_Phred', Number=>1, Type=>'Float', Description=>'CADD 1.3 Phred Score (Combined Annotation Dependent Depletion) is a score that is based on SVM on multiple other scores. It assigns a score to each possible mutation in the human genome, therefore can evaluate non-coding variants as well as coding ones - version 1.3 - Annovar 20170123'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype cadd13
            $VERBOSE1 and printerr("\n**********\n* Resuming WGSCADD annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype cadd13 -buildver $USEDBUILD -otherinfo -indexfilter_threshold 1 -outfile $OUTFILE_WGSCADD $ANNOINFILE $WGSCADDPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of WGSCADD annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of WGSCADD annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_WGSCADD, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";
        
        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_WGSCADD>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(cadd13)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($cadd, $caddphred) = split(/,/, $2);                
		$rH_data->{$key}->{WGSCADD} = $cadd;
                $rH_data->{$key}->{WGSCADD_Phred} = $caddphred;
                
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_WGSCADD);
    }

    # if not excluded, run ANNOVAR -filter to output WGSFATHMM annotations
    if (doIt(-elem => "wgsfathmm")) {
        my $annovarOutput = "$OUTFILE_WGSFATHMM.$USEDBUILD\_fathmm_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'WGSFATHMM_coding', Number=>1, Type=>'Float', Description=>'whole-genome FATHMM_coding Categorical Prediction score, Very similar to SIFT and Polyphen. If a score is smaller than -1.5 the corresponding NS is predicted as D(AMAGING); otherwise it is predicted as T(OLERATED) - version_20160315'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'WGSFATHMM_noncoding', Number=>1, Type=>'Float', Description=>'whole-genome FATHMM_noncoding score - version_20160315'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype fathmm
            $VERBOSE1 and printerr("\n**********\n* Resuming WGSFATHMM annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype fathmm -buildver $USEDBUILD -otherinfo -indexfilter_threshold 1 -outfile $OUTFILE_WGSFATHMM $ANNOINFILE $WGSFATHMMPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of WGSFATHMM annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of WGSFATHMM annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_WGSFATHMM, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_WGSFATHMM>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(fathmm)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($fathmm_noncoding, $fathmm_coding) = split(/,/, $2);
                $rH_data->{$key}->{WGSFATHMM_coding} = $fathmm_coding;
                $rH_data->{$key}->{WGSFATHMM_noncoding} = $fathmm_noncoding;

            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_WGSFATHMM);
    }

    # if not excluded, run ANNOVAR -filter to output HRCR1 annotations
    if (doIt(-elem => "hrcr1")) {
        my $annovarOutput = "$OUTFILE_HRCR1.$USEDBUILD\_hrcr1_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'HRCR1_AF', Number=>1, Type=>'Float', Description=>'Allele frequency in hrcr1 database (40 million variants from 32K samples in haplotype reference consortium) - version_20151203'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'HRCR1_AC', Number=>1, Type=>'Float', Description=>'Allele count in hrcr1 database (40 million variants from 32K samples in haplotype reference consortium) - version_20151203'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'HRCR1_AN', Number=>1, Type=>'Float', Description=>'Allele number in hrcr1 database (40 million variants from 32K samples in haplotype reference consortium) - version_20151203'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'HRCR1_nonKG_AF', Number=>1, Type=>'Float', Description=>'Non KG Allele frequency in hrcr1 database (40 million variants from 32K samples in haplotype reference consortium) - version_20151203'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'HRCR1_nonKG_AC', Number=>1, Type=>'Float', Description=>'Non KG Allele count in hrcr1 database (40 million variants from 32K samples in haplotype reference consortium) - version_20151203'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'HRCR1_nonKG_AN', Number=>1, Type=>'Float', Description=>'Non KG Allele number in hrcr1 database (40 million variants from 32K samples in haplotype reference consortium - version_20151203'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype hrcr1
            $VERBOSE1 and printerr("\n**********\n* Resuming HRCR1 annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype hrcr1 -buildver $USEDBUILD -otherinfo -indexfilter_threshold 1 -outfile $OUTFILE_HRCR1 $ANNOINFILE $HRCR1PATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of HRCR1 annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of HRCR1 annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_HRC, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
	$VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_HRC>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(hrcr1)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($hrc_af, $hrc_ac, $hrc_an, $hrc_non1kg_af, $hrc_non1kg_ac, $hrc_non1kg_an) = split(/,/, $2);
                $rH_data->{$key}->{HRCR1_AF}       = $hrc_af;
                $rH_data->{$key}->{HRCR1_AC}       = $hrc_ac;
                $rH_data->{$key}->{HRCR1_AN}       = $hrc_an;
                $rH_data->{$key}->{HRCR1_nonKG_AF} = $hrc_non1kg_af;
                $rH_data->{$key}->{HRCR1_nonKG_AC} = $hrc_non1kg_ac;
                $rH_data->{$key}->{HRCR1_nonKG_AN} = $hrc_non1kg_an;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_HRC);
    }

    # if not excluded, run ANNOVAR -filter to output DBSCSNV11 annotations
    if (doIt(-elem => "dbscsnv11")) {
        my $annovarOutput = "$OUTFILE_DBSCSNV11.$USEDBUILD\_dbscsnv11_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbscSNV_ADA_SCORE', Number=>1, Type=>'Float', Description=>'It provides splice site effect prediction by AdaBoost and Random Forest from dbscSNV version 1.1 - 20151218'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'dbscSNV_RF_SCORE', Number=>1, Type=>'Float', Description=>'It provides splice site effect prediction by AdaBoost and Random Forest from dbscSNV version 1.1 - 20151218'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype dbscsnv11
            $VERBOSE1 and printerr("\n**********\n* Resuming DBSCSNV11 annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype dbscsnv11 -buildver $USEDBUILD -otherinfo -indexfilter_threshold 1 -outfile $OUTFILE_DBSCSNV11 $ANNOINFILE $DBSCSNVPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of DBSCSNV11 annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of DBSCSNV11 annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_DBSCSNV11, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_DBSCSNV11>) {         # find value of Lookup in ANNOVAR output, and annotate variant in %results
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(dbscsnv11)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my ($adascore, $rfscore) = split(/,/, $2); #dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE
                $rH_data->{$key}->{dbscSNV_ADA_SCORE} = $adascore;
                $rH_data->{$key}->{dbscSNV_RF_SCORE} = $rfscore;

            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSCSNV11);
    }

    # if not excluded, run ANNOVAR -filter to output EXAC_ALL_nonpsych annotations
    if (doIt(-elem => "exac_nonpsych")) {
        my $annovarOutput = "$OUTFILE_ExAC_ALL_nonpsych.$USEDBUILD\_exac03nonpsych_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_ALL_nonpsych', Number=>1, Type=>'Float', Description=>'ExAC on non-Psychiatric disease samples allele frequency - version_20151203'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype exac_all_nonpsych
            $VERBOSE1 and printerr("\n**********\n* Resuming ExAC_all_nonpsych annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype exac03nonpsych -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ExAC_ALL_nonpsych $ANNOINFILE $EXACPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ExAC_all_nonpsych annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of EXAC_all_nonpsych annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ExAC_nonpsych, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ExAC_nonpsych>) {
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(exac03nonpsych)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ExAC_ALL_nonpsych} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ExAC_nonpsych);
    }

    # if not excluded, run ANNOVAR -filter to output EXAC_ALL_nontcga annotations
    if (doIt(-elem => "exac_nontcga")) {
        my $annovarOutput = "$OUTFILE_ExAC_ALL_nontcga.$USEDBUILD\_exac03nontcga_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_ALL_nontcga', Number=>1, Type=>'Float', Description=>'ExAC on non-TCGA samples allele frequency - version_20151203'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype exac_all_nontcga
            $VERBOSE1 and printerr("\n**********\n* Resuming ExAC_all_nontcga annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype exac03nontcga -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ExAC_ALL_nontcga $ANNOINFILE $EXACPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ExAC_all_nontcga annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of EXAC_all_nontcga annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ExAC_nontcga, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ExAC_nontcga>) {
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(exac03nontcga)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ExAC_ALL_nontcga} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ExAC_nontcga);
    }

    # if not excluded, run ANNOVAR -filter to output kaviar annotations
    if (doIt(-elem => "kaviar")) {
        my $annovarOutput = "$OUTFILE_KAVIAR.$USEDBUILD\_kaviar_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KAVIAR_AF', Number=>1, Type=>'Float', Description=>'Kaviar Allele Frequency based on 170 million Known Variants from 13K genomes and 64K exomes in 34 projects - updated_20151203'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KAVIAR_AC', Number=>1, Type=>'Float', Description=>'Kaviar Allele Count based on 170 million Known Variants from 13K genomes and 64K exomes in 34 projects - updated_20151203'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'KAVIAR_AN', Number=>1, Type=>'Float', Description=>'Kaviar Allele Number based on 170 million Known Variants from 13K genomes and 64K exomes in 34 projects - updated_20151203'});

        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype kaviar
            $VERBOSE1 and printerr("\n**********\n* Resuming kaviar annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype kaviar -buildver $USEDBUILD -otherinfo -indexfilter_threshold 1 -outfile $OUTFILE_KAVIAR $ANNOINFILE $KAVIARPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of kaviar annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of kaviar annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_KAVIAR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_KAVIAR>) {
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(kaviar)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                 my ($aF, $aC, $aN) = split(/,/, $2);
                $rH_data->{$key}->{KAVIAR_AF} = $aF;
                $rH_data->{$key}->{KAVIAR_AC} = $aC;
                $rH_data->{$key}->{KAVIAR_AN} = $aN;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_KAVIAR);
    }

    # if not excluded, run ANNOVAR -filter to output WHOLE GENOME EIGEN annotations
    if (doIt(-elem => "wgseigen")) {
        my $annovarOutput = "$OUTFILE_WGSEIGEN.$USEDBUILD\_wgseigen_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'WGSEIGEN', Number=>1, Type=>'Float', Description=>'Whole-genome Eigen scores, a spectral approach integrating functional genomic annotations for coding and noncoding variants - version_20160330'});

        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype wgseigen
            $VERBOSE1 and printerr("\n**********\n* Resuming WHOLE GENOME EIGEN annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype wgseigen -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_WGSEIGEN $ANNOINFILE $WGSEIGENPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of WGS EIGEN annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of WGS EIGEN annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_WGSEIGEN, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_WGSEIGEN>) {
            chomp $line;    # $1      $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(wgseigen)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{WGSEIGEN} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_WGSEIGEN);
    }


    # if not excluded, run ANNOVAR -filter to output NCI60 annotations
    if (doIt(-elem => "nci60")) {
        my $annovarOutput = "$OUTFILE_NCI60.$USEDBUILD\_nci60_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'NCI60', Number=>'A', Type=>'String', Description=>'NCI-60 human tumor cell line panel exome sequencing allele frequency data (ANNOVAR - nci60)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype nci60
            $VERBOSE1 and printerr("\n**********\n* Resuming NCI60 annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype nci60 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_NCI60 $ANNOINFILE $CANCERPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of NCI60 annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of NCI60 annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_NCI60, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_NCI60>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1    $2     $3    $4        $5     $6       $7         $8
            if ($line =~ /(nci60)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{NCI60} = $2;
                $rH_data->{$key}->{NCI60} =~ s/\=/\:/g;
                $rH_data->{$key}->{NCI60} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{NCI60} =~ s/ /_/g;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_NCI60);
    }

    # if not excluded, run ANNOVAR -filter to output GoNL Frequencies
    if (doIt(-elem => "gonl")) {
        my $annovarOutput = "$OUTFILE_GONL.$USEDBUILD\_vcf_dropped";
        my $vcfdbfile = "gonl.snp_indels.r5.vcf";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GoNL', Number=>'A', Type=>'String', Description=>'GoNL project : frequencies from 499 unrelated individuals - release 5'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype vcf
            $VERBOSE1 and printerr("\n**********\n* Resuming GoNL annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype vcf -vcfdbfile $vcfdbfile -infoasscore -buildver $USEDBUILD -outfile $OUTFILE_GONL $ANNOINFILE $FREQUENCIESPATH/GoNL/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GoNL annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GoNL annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GONL, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GONL>) {
            chomp $line;  #  $1   $2     $3   $4         $5     $6       $7           $8
            if ($line =~ /(vcf)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                ($rH_data->{$key}->{GoNL} = $2) =~ s/.+AF=([^;]+).+/$1/;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GONL);
    }

    # if not excluded, run ANNOVAR -filter to output the frequencies from the Wellderly project
    if (doIt(-elem => "wellderly")) {
        my $annovarOutput = "$OUTFILE_WELLDERLY.$USEDBUILD\_vcf_dropped";
        my $vcfdbfile = "SWGR_v1.0_snvs_small_indels_shuffled_ALL.vcf";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'WELLDERLY_AF', Number=>'A', Type=>'String', Description=>'Allele Frequencies over 454 unrelated Scripps Wellderly Study participants with European ancestry - version SWGR_v1.0'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'WELLDERLY_AN', Number=>'A', Type=>'String', Description=>'Total Allele Counts over 454 unrelated Scripps Wellderly Study participants with European ancestry - version SWGR_v1.0'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype vcf
            $VERBOSE1 and printerr("\n**********\n* Resuming Wellderly annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype vcf -vcfdbfile $vcfdbfile -infoasscore -buildver $USEDBUILD -outfile $OUTFILE_WELLDERLY $ANNOINFILE $FREQUENCIESPATH/Wellderly/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of Wellderly annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of Wellderly annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_WELLDERLY, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_WELLDERLY>) {
            chomp $line; #  $1    $2     $3   $4         $5     $6        $7          $8
            if ($line =~ /(vcf)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $annot = $2;
                ($rH_data->{$key}->{WELLDERLY_AF} = $annot) =~ s/.*AF=([^;]+).*/$1/;
                ($rH_data->{$key}->{WELLDERLY_AN} = $annot) =~ s/.*AN=([^;]+).*/$1/;
#                my ($AF, $AN) = split(/\|/, $2);
#                ($rH_data->{$key}->{WELLDERLY_AF} = $AF) =~ s/^AF://;
#                ($rH_data->{$key}->{WELLDERLY_AN} = $AN) =~ s/^AN://;
            } else {
print "after if (fails): line = \"$line\"\n";
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_WELLDERLY);
    }

########################################################################################################################################################
# DBNSFP start
    # if not excluded, run ANNOVAR -filter to output DBNSFP annotations
    if (doIt(-elem => "dbnsfp") || doIt(-elem => "sift") || doIt(-elem => "pp2") || doIt(-elem => "lrt") || doIt(-elem => "mt")
            || doIt(-elem => "ma") || doIt(-elem => "fathmm") || doIt(-elem => "provean") || doIt(-elem => "vest") || doIt(-elem => "msvm") || doIt(-elem => "mlr") || doIt(-elem => "mcap") || doIt(-elem => "cadd") || doIt(-elem => "dann") || doIt(-elem => "fathmmkl") || doIt(-elem => "eigen") || doIt(-elem => "genocanyon") || doIt(-elem => "intgr") || doIt(-elem => "gerp") || doIt(-elem => "phylop") || doIt(-elem => "phcons") || doIt(-elem => "siphy") || doIt(-elem => "interpro") || doIt(-elem => "gtex") ) {
#print STDERR "doIt = yes\n";
        my $annovarOutput = "$OUTFILE_DBNSFP.$USEDBUILD\_dbnsfp33a_dropped";
        (doIt(-elem => "sift"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'SIFT_score', Number=>'A', Type=>'Float', Description=>'whole-exome SIFT scores with missing values imputed (dbnsfp version 3.3a - 20170221)'});
	(doIt(-elem => "sift"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'SIFT_rank', Number=>'A', Type=>'Float', Description=>'whole-exome SIFT rank with missing values imputed (dbnsfp version 3.3a - 20170221)'});
        (doIt(-elem => "sift"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'SIFT_pred', Number=>'A', Type=>'String', Description=>'whole-exome SIFT categorical prediction with missing values imputed (dbnsfp version 3.3a - 20170221)'});
        (doIt(-elem => "pp2"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'Polyphen2_HDIV_score', Number=>'A', Type=>'Float', Description=>'whole-exome PolyPhen 2 scores built on HumanDiv database (for complex phenotypes) (dbnsfp version 3.3a - 20170221)'});
	(doIt(-elem => "pp2"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'Polyphen2_HDIV_rank', Number=>'A', Type=>'Float', Description=>'whole-exome PolyPhen 2 rank built on HumanDiv database (for complex phenotypes) (dbnsfp version 3.3a - 20170221)'});
        (doIt(-elem => "pp2"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'Polyphen2_HDIV_pred', Number=>'A', Type=>'String', Description=>'whole-exome PolyPhen 2 categorical prediction built on HumanDiv database (for complex phenotypes) (dbnsfp version 3.3a - 20170221)'});
        (doIt(-elem => "pp2"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'Polyphen2_HVAR_score', Number=>'A', Type=>'Float', Description=>'whole-exome PolyPhen 2 scores built on HumanVar database (for Mendelian phenotypes) (dbnsfp version 3.3a - 20170221)'});
        (doIt(-elem => "pp2"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'Polyphen2_HVAR_rank', Number=>'A', Type=>'Float', Description=>'whole-exome PolyPhen 2 rank built on HumanVar database (for Mendelian phenotypes) (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "pp2"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'Polyphen2_HVAR_pred', Number=>'A', Type=>'String', Description=>'whole-exome PolyPhen 2 categorical prediction built on HumanVar database (for Mendelian phenotypes) (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "lrt"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'LRT_score', Number=>'A', Type=>'Float', Description=>'whole-exome LRT scores inccluding raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "lrt"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'LRT_rank', Number=>'A', Type=>'Float', Description=>'whole-exome LRT rank inccluding raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "lrt"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'LRT_pred', Number=>'A', Type=>'String', Description=>'whole-exome LRT scores inccluding categorical prediction score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mt"))       and $rO_vcf->add_header_line({key=>'INFO', ID=>'MutationTaster_score', Number=>'A', Type=>'Float', Description=>'whole-exome MutationTaster scores inccluding raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mt"))       and $rO_vcf->add_header_line({key=>'INFO', ID=>'MutationTaster_rank', Number=>'A', Type=>'Float', Description=>'whole-exome MutationTaster rank inccluding raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mt"))       and $rO_vcf->add_header_line({key=>'INFO', ID=>'MutationTaster_pred', Number=>'A', Type=>'String', Description=>'whole-exome MutationTaster categorical prediction inccluding raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "ma"))       and $rO_vcf->add_header_line({key=>'INFO', ID=>'MutationAssessor_score', Number=>'A', Type=>'Float', Description=>'whole-exome MutationAssessor scores including raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "ma"))       and $rO_vcf->add_header_line({key=>'INFO', ID=>'MutationAssessor_rank', Number=>'A', Type=>'Float', Description=>'whole-exome MutationAssessor rank including raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "ma"))       and $rO_vcf->add_header_line({key=>'INFO', ID=>'MutationAssessor_pred', Number=>'A', Type=>'String', Description=>'whole-exome MutationAssessor categorical prediction inccluding raw score and categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "fathmm"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'FATHMM_score', Number=>'A', Type=>'Float', Description=>'whole-exome FATHMM scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "fathmm"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'FATHMM_rank', Number=>'A', Type=>'Float', Description=>'whole-exome FATHMM rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "fathmm"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'FATHMM_pred', Number=>'A', Type=>'String', Description=>'whole-exome FATHMM categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "provean"))  and $rO_vcf->add_header_line({key=>'INFO', ID=>'PROVEAN_score', Number=>'A', Type=>'Float', Description=>'whole-exome PROVEAN scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "provean"))  and $rO_vcf->add_header_line({key=>'INFO', ID=>'PROVEAN_rank', Number=>'A', Type=>'Float', Description=>'whole-exome PROVEAN rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "provean"))  and $rO_vcf->add_header_line({key=>'INFO', ID=>'PROVEAN_pred', Number=>'A', Type=>'Float', Description=>'whole-exome PROVEAN categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "vest"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'VEST3_score', Number=>'A', Type=>'Float', Description=>'whole-exome VEST scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "vest"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'VEST3_rank', Number=>'A', Type=>'Float', Description=>'whole-exome VEST rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "msvm"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'MetaSVM_score', Number=>'A', Type=>'Float', Description=>'whole-exome MetaSVM scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "msvm"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'MetaSVM_rank', Number=>'A', Type=>'Float', Description=>'whole-exome MetaSVM rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "msvm"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'MetaSVM_pred', Number=>'A', Type=>'String', Description=>'whole-exome MetaSVM categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mlr"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'MetaLR_score', Number=>'A', Type=>'Float', Description=>'whole-exome MetaLR scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mlr"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'MetaLR_rank', Number=>'A', Type=>'Float', Description=>'whole-exome MetaLR rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mlr"))      and $rO_vcf->add_header_line({key=>'INFO', ID=>'MetaLR_pred', Number=>'A', Type=>'String', Description=>'whole-exome MetaLR categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mcap"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'M_CAP_score', Number=>'A', Type=>'Float', Description=>'whole-exome M-CAP scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mcap"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'M_CAP_rank', Number=>'A', Type=>'Float', Description=>'whole-exome M-CAP rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "mcap"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'M_CAP_pred', Number=>'A', Type=>'String', Description=>'whole-exome M-CAP categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "cadd"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'CADD_raw', Number=>'A', Type=>'Float', Description=>'whole-exome CADD raw scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "cadd"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'CADD_rank', Number=>'A', Type=>'Float', Description=>'whole-exome CADD rank scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "cadd"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'CADD_phred', Number=>'A', Type=>'Float', Description=>'whole-exome CADD phred scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "dann"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'DANN_score', Number=>'A', Type=>'Float', Description=>'whole-exome DANN scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "dann"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'DANN_rank', Number=>'A', Type=>'Float', Description=>'whole-exome DANN rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "fathmmkl")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'FATHMMMKL_coding_score', Number=>'A', Type=>'Float', Description=>'whole-exome FATHMM scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "fathmmkl")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'FATHMMMKL_coding_rank', Number=>'A', Type=>'Float', Description=>'whole-exome FATHMM rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "fathmmkl")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'FATHMMMKL_coding_pred', Number=>'A', Type=>'String', Description=>'whole-exome FATHMM categorical prediction (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "eigen"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'EIGEN_coding_or_noncoding', Number=>'A', Type=>'String', Description=>'whole-exome EIGEN coding or non-coding (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "eigen"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'EIGEN_raw', Number=>'A', Type=>'Float', Description=>'whole-exome EIGEN raw score (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "eigen"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'EIGEN_PC_raw', Number=>'A', Type=>'Float', Description=>'whole-exome EIGEN PC raw score (dbnsfp version 3.3a - 20170221)'});
	      (doIt(-elem => "genocanyon"))  and $rO_vcf->add_header_line({key=>'INFO', ID=>'GenoCanyon_score', Number=>'A', Type=>'Float', Description=>'whole-exome GenoCanyon scores (dbnsfp version 3.3a - 20170221)'});
	      (doIt(-elem => "genocanyon"))  and $rO_vcf->add_header_line({key=>'INFO', ID=>'GenoCanyon_rank', Number=>'A', Type=>'Float', Description=>'whole-exome GenoCanyon rank (dbnsfp version 3.3a - 20170221)'});
	       (doIt(-elem => "intgr"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'INTEGRATED_fCons_score', Number=>'A', Type=>'Float', Description=>'whole-exome Integrated fitCons score (dbnsfp version 3.3a - 20170221)'});   
		(doIt(-elem => "intgr"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'INTEGRATED_fCons_rank', Number=>'A', Type=>'Float', Description=>'whole-exome Integrated fitCons rank (dbnsfp version 3.3a - 20170221)'});
	       (doIt(-elem => "intgr"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'INTEGRATED_conf_value', Number=>'A', Type=>'Float', Description=>'whole-exome Integrated confidence value (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "gerp"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'GERP_RS_dbnsfp', Number=>'A', Type=>'Float', Description=>'whole-exome GERP++ scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "gerp"))     and $rO_vcf->add_header_line({key=>'INFO', ID=>'GERP_rank_dbnsfp', Number=>'A', Type=>'Float', Description=>'whole-exome GERP++ rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phylop"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phyloP100way_vertebrate_score', Number=>'A', Type=>'Float', Description=>'whole-exome PhyloP scores based on 100-way alignment vertebrate subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phylop"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phyloP100way_vertebrate_rank', Number=>'A', Type=>'Float', Description=>'whole-exome PhyloP rank based on 100-way alignment vertebrate subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phylop"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phyloP20way_mammalian_score', Number=>'A', Type=>'Float', Description=>'whole-exome PhyloP scores based on 20-way alignment mammalian subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phylop"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phyloP20way_mammalian_rank', Number=>'A', Type=>'Float', Description=>'whole-exome PhyloP rank based on 20-way alignment mammalian subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phcons"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phastCons100way_vertebrate_score', Number=>'A', Type=>'Float', Description=>'whole-exome PhastCons scores based on 100-way alignment vertebrate subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phcons"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phastCons100way_vertebrate_rank', Number=>'A', Type=>'Float', Description=>'whole-exome PhastCons rank based on 100-way alignment vertebrate subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phcons"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phastCons20way_mammalian_score', Number=>'A', Type=>'Float', Description=>'whole-exome PhastCons scores based on 100-way alignment mammalian subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "phcons"))   and $rO_vcf->add_header_line({key=>'INFO', ID=>'phastCons20way_mammalian_rank', Number=>'A', Type=>'Float', Description=>'whole-exome PhastCons rank based on 100-way alignment mammalian subset (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "siphy"))    and $rO_vcf->add_header_line({key=>'INFO', ID=>'SiPhy_29way_logOdds_score', Number=>'A', Type=>'Float', Description=>'whole-exome SiPhy scores (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "siphy"))    and $rO_vcf->add_header_line({key=>'INFO', ID=>'SiPhy_29way_logOdds_rank', Number=>'A', Type=>'Float', Description=>'whole-exome SiPhy rank (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "interpro")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'Interpro', Number=>'A', Type=>'String', Description=>'whole-exome Interpro domains (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "gtex")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'GTEx_V6_gene', Number=>'A', Type=>'String', Description=>'whole-exome Genotype-Tissue Expression (GTex) V6 Gene (dbnsfp version 3.3a - 20170221)'});
		(doIt(-elem => "gtex")) and $rO_vcf->add_header_line({key=>'INFO', ID=>'GTEx_V6_tissue', Number=>'A', Type=>'String', Description=>'whole-exome Genotype-Tissue Expression (GTex) V6 Tissue (dbnsfp version 3.3a - 20170221)'});

		if ($RESUMEFOLDER && -e $annovarOutput) {
		    # resuming ANNOVAR -filter -dbtype dbnsfp
		    $VERBOSE1 and printerr("\n**********\n* Resuming DBNSFP annotations using already existing file $annovarOutput\n**********\n\n");
		} else {
		    my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype dbnsfp33a -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBNSFP $ANNOINFILE $DBNSFPPATH/ $VERBOSE1";
		    $VERBOSE1 and printerr("\n**********\n* Starting fetching of DBNSFP annotations with system command:\n* <$systemCommand>\n**********\n\n");
		    my $error = system ($systemCommand);
		    if ($error) {
			printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
		    } elsif ($VERBOSE2) {
			my $exitSignal = ($? >> 8);
			my $message = sprintf "Command exited with value %d", $exitSignal;
			printerr("\n**********\n* $message\n**********\n\n");
		    }
		    $VERBOSE1 and printerr("\n**********\n* Fetching of DBNSFP annotations scores done...\n**********\n\n");
		}

		# open ANNOVAR output file for processing
		open(ANNOVAR_DBNSFP, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

		# annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
		$VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
		while (my $line = <ANNOVAR_DBNSFP>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
		    chomp $line;    #  $1       $2     $3    $4        $5     $6       $7          $8
	#              dbnsfp33a	0.027,D,0.146,B,0.101,B,0.000,D,0.000,P,2.075,M,-0.88,T,-2.69,D,0.369,2.294,18.13,0.972,0.978,D,-0.592,T,0.282,T,0.732,0,5.2,0.871,0.883,0.379,0.361,19.102	1	100358103	100358103	C	T	1	100358103	C	T
		    if ($line =~ /(dbnsfp33a)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
			my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
			my ($SIFT_score, $SIFT_rank, $SIFT_pred, 
			    $Polyphen2_HDIV_score, $Polyphen2_HDIV_rank, $Polyphen2_HDIV_pred, $Polyphen2_HVAR_score, $Polyphen2_HVAR_rank, $Polyphen2_HVAR_pred, 
			    $LRT_score, $LRT_rank, $LRT_pred, 
			    $MutationTaster_score, $MutationTaster_rank, $MutationTaster_pred, 
			    $MutationAssessor_score, $MutationAssessor_rank, $MutationAssessor_pred, 
			    $FATHMM_score, $FATHMM_rank, $FATHMM_pred, 
			    $PROVEAN_score, $PROVEAN_rank, $PROVEAN_pred, 
			    $VEST3_score, $VEST3_rank, 
			    $MetaSVM_score, $MetaSVM_rank, $MetaSVM_pred, 
			    $MetaLR_score, $MetaLR_rank, $MetaLR_pred, 
			    $M_CAP_score, $M_CAP_rank, $M_CAP_pred, 
			    $CADD_raw, $CADD_rank, $CADD_phred, 
			    $DANN_score, $DANN_rank, 
			    $FATHMMMKL_coding_score, $FATHMMMKL_coding_rank, $FATHMMMKL_coding_pred, 
			    $EIGEN_coding_or_noncoding, $EIGEN_raw, $EIGEN_PC_raw, 
			    $GenoCanyon_score, $GenoCanyon_rank,
			    $INTEGRATED_fCons_score, $INTEGRATED_fCons_rank, $INTEGRATED_conf_value, 
			    $GERP_RS_dbnsfp, $GERP_rank_dbnsfp, 
			    $phyloP100way_vertebrate_score, $phyloP100way_vertebrate_rank, $phyloP20way_mammalian_score, $phyloP20way_mammalian_rank, 
			    $phastCons100way_vertebrate_score, $phastCons100way_vertebrate_rank, $phastCons20way_mammalian_score, $phastCons20way_mammalian_rank, 
			    $SiPhy_29way_logOdds_score, $SiPhy_29way_logOdds_rank, 
			    $Interpro, 
			    $GTEx_V6_gene, $GTEx_V6_tissue) = split /,/, $2;

			if ( $Interpro ) { 
				$Interpro =~ s/ /_/g;
				$Interpro =~ s/;/|/g;
				$Interpro =~ s/__/_/g;
	#                        print  $Interpro . "\n";
			}
	#                 print $key ."\n";

			(doIt(-elem => "sift"))   and $rH_data->{$key}->{SIFT_score}                = ($SIFT_score eq ".")                ? undef : $SIFT_score;
			(doIt(-elem => "sift"))   and $rH_data->{$key}->{SIFT_rank}                 = ($SIFT_rank  eq ".")                ? undef : $SIFT_rank;
			(doIt(-elem => "sift"))   and $rH_data->{$key}->{SIFT_pred}                 = ($SIFT_pred eq ".")                 ? undef : $SIFT_pred;
			(doIt(-elem => "pp2"))    and $rH_data->{$key}->{Polyphen2_HDIV_score}      = ($Polyphen2_HDIV_score eq ".")      ? undef : $Polyphen2_HDIV_score;
                (doIt(-elem => "pp2"))    and $rH_data->{$key}->{Polyphen2_HDIV_rank}       = ($Polyphen2_HDIV_rank eq ".")       ? undef : $Polyphen2_HDIV_rank;
                (doIt(-elem => "pp2"))    and $rH_data->{$key}->{Polyphen2_HDIV_pred}       = ($Polyphen2_HDIV_pred eq ".")       ? undef : $Polyphen2_HDIV_pred;
                (doIt(-elem => "pp2"))    and $rH_data->{$key}->{Polyphen2_HVAR_score}      = ($Polyphen2_HVAR_score eq ".")      ? undef : $Polyphen2_HVAR_score;
                (doIt(-elem => "pp2"))    and $rH_data->{$key}->{Polyphen2_HVAR_rank}       = ($Polyphen2_HVAR_rank eq ".")       ? undef : $Polyphen2_HVAR_rank;
                (doIt(-elem => "pp2"))    and $rH_data->{$key}->{Polyphen2_HVAR_pred}       = ($Polyphen2_HVAR_pred eq ".")       ? undef : $Polyphen2_HVAR_pred;
                (doIt(-elem => "lrt"))    and $rH_data->{$key}->{LRT_score}                 = ($LRT_score eq ".")                 ? undef : $LRT_score;
                (doIt(-elem => "lrt"))    and $rH_data->{$key}->{LRT_rank}                  = ($LRT_rank eq ".")                  ? undef : $LRT_rank;
                (doIt(-elem => "lrt"))    and $rH_data->{$key}->{LRT_pred}                  = ($LRT_pred eq ".")                  ? undef : $LRT_pred;
                (doIt(-elem => "mt"))     and $rH_data->{$key}->{MutationTaster_score}      = ($MutationTaster_score eq ".")      ? undef : $MutationTaster_score;
                (doIt(-elem => "mt"))     and $rH_data->{$key}->{MutationTaster_rank}       = ($MutationTaster_rank eq ".")       ? undef : $MutationTaster_rank;
                (doIt(-elem => "mt"))     and $rH_data->{$key}->{MutationTaster_pred}       = ($MutationTaster_pred eq ".")       ? undef : $MutationTaster_pred;
                (doIt(-elem => "ma"))     and $rH_data->{$key}->{MutationAssessor_score}    = ($MutationAssessor_score eq ".")    ? undef : $MutationAssessor_score;
                (doIt(-elem => "ma"))     and $rH_data->{$key}->{MutationAssessor_rank}     = ($MutationAssessor_rank eq ".")     ? undef : $MutationAssessor_rank;
                (doIt(-elem => "ma"))     and $rH_data->{$key}->{MutationAssessor_pred}     = ($MutationAssessor_pred eq ".")     ? undef : $MutationAssessor_pred;
                (doIt(-elem => "fathmm")) and $rH_data->{$key}->{FATHMM_score}              = ($FATHMM_score eq ".")              ? undef : $FATHMM_score;
                (doIt(-elem => "fathmm")) and $rH_data->{$key}->{FATHMM_rank}               = ($FATHMM_rank eq ".")               ? undef : $FATHMM_rank;
                (doIt(-elem => "fathmm")) and $rH_data->{$key}->{FATHMM_pred}               = ($FATHMM_pred eq ".")               ? undef : $FATHMM_pred;
		(doIt(-elem => "provean")) and $rH_data->{$key}->{PROVEAN_score}            = ($PROVEAN_score eq ".")             ? undef : $PROVEAN_score;
                (doIt(-elem => "provean")) and $rH_data->{$key}->{PROVEAN_rank}             = ($PROVEAN_rank eq ".")              ? undef : $PROVEAN_rank;
		(doIt(-elem => "provean")) and $rH_data->{$key}->{PROVEAN_pred}             = ($PROVEAN_pred eq ".")              ? undef : $PROVEAN_pred;
                (doIt(-elem => "vest"))   and $rH_data->{$key}->{VEST3_score}               = ($VEST3_score eq ".")               ? undef : $VEST3_score;
                (doIt(-elem => "vest"))   and $rH_data->{$key}->{VEST3_rank}                = ($VEST3_rank eq ".")                ? undef : $VEST3_rank;
                (doIt(-elem => "msvm")) and $rH_data->{$key}->{MetaSVM_score}               = ($MetaSVM_score eq ".")             ? undef : $MetaSVM_score;
                (doIt(-elem => "msvm")) and $rH_data->{$key}->{MetaSVM_rank}                = ($MetaSVM_rank eq ".")              ? undef : $MetaSVM_rank;
                (doIt(-elem => "msvm")) and $rH_data->{$key}->{MetaSVM_pred}                = ($MetaSVM_pred eq ".")              ? undef : $MetaSVM_pred;
                (doIt(-elem => "mlr")) and $rH_data->{$key}->{MetaLR_score}                 = ($MetaLR_score eq ".")              ? undef : $MetaLR_score;
                (doIt(-elem => "mlr")) and $rH_data->{$key}->{MetaLR_rank}                  = ($MetaLR_rank eq ".")               ? undef : $MetaLR_rank;
                (doIt(-elem => "mlr")) and $rH_data->{$key}->{MetaLR_pred}                  = ($MetaLR_pred eq ".")               ? undef : $MetaLR_pred;
                (doIt(-elem => "mcap")) and $rH_data->{$key}->{M_CAP_score}                 = ($M_CAP_score eq ".")               ? undef : $M_CAP_score;
                (doIt(-elem => "mcap")) and $rH_data->{$key}->{M_CAP_rank}                  = ($M_CAP_rank eq ".")                ? undef : $M_CAP_rank;
                (doIt(-elem => "mcap")) and $rH_data->{$key}->{M_CAP_pred}                  = ($M_CAP_pred eq ".")                ? undef : $M_CAP_pred;
                (doIt(-elem => "cadd"))   and $rH_data->{$key}->{CADD_raw}                  = ($CADD_raw eq ".")                  ? undef : $CADD_raw;
                (doIt(-elem => "cadd"))   and $rH_data->{$key}->{CADD_rank}                 = ($CADD_rank eq ".")                 ? undef : $CADD_rank;
                (doIt(-elem => "cadd"))   and $rH_data->{$key}->{CADD_phred}                = ($CADD_phred eq ".")                ? undef : $CADD_phred;
		(doIt(-elem => "dann")) and $rH_data->{$key}->{DANN_score}                  = ($DANN_score eq ".")                ? undef : $DANN_score; 
                (doIt(-elem => "dann")) and $rH_data->{$key}->{DANN_rank}                   = ($DANN_rank eq ".")                 ? undef : $DANN_rank;
		(doIt(-elem => "fathmmkl")) and $rH_data->{$key}->{FATHMMMKL_coding_score}  = ($FATHMMMKL_coding_score eq ".")    ? undef : $FATHMMMKL_coding_score;
                (doIt(-elem => "fathmmkl")) and $rH_data->{$key}->{FATHMMMKL_coding_rank}   = ($FATHMMMKL_coding_rank eq ".")     ? undef : $FATHMMMKL_coding_rank;
		(doIt(-elem => "fathmmkl")) and $rH_data->{$key}->{FATHMMMKL_coding_pred}   = ($FATHMMMKL_coding_pred eq ".")     ? undef : $FATHMMMKL_coding_pred;
                (doIt(-elem => "eigen")) and $rH_data->{$key}->{EIGEN_coding_or_noncoding}  = ($EIGEN_coding_or_noncoding eq ".") ? undef : $EIGEN_coding_or_noncoding;
                (doIt(-elem => "eigen")) and $rH_data->{$key}->{EIGEN_raw}                  = ($EIGEN_raw eq ".")                 ? undef : $EIGEN_raw;
                (doIt(-elem => "eigen")) and $rH_data->{$key}->{EIGEN_PC_raw}               = ($EIGEN_PC_raw eq ".")              ? undef : $EIGEN_PC_raw;
                (doIt(-elem => "genocanyon")) and $rH_data->{$key}->{GenoCanyon_score}      = ($GenoCanyon_score eq ".")          ? undef : $GenoCanyon_score;
                (doIt(-elem => "genocanyon")) and $rH_data->{$key}->{GenoCanyon_rank}       = ($GenoCanyon_rank eq ".")           ? undef : $GenoCanyon_rank;
		(doIt(-elem => "intgr")) and $rH_data->{$key}->{INTEGRATED_fCons_score}     = ($INTEGRATED_fCons_score eq ".")    ? undef : $INTEGRATED_fCons_score;
                (doIt(-elem => "intgr")) and $rH_data->{$key}->{INTEGRATED_fCons_rank}      = ($INTEGRATED_fCons_rank eq ".")     ? undef : $INTEGRATED_fCons_rank;
		(doIt(-elem => "intgr")) and $rH_data->{$key}->{INTEGRATED_conf_value}      = ($INTEGRATED_conf_value eq ".")     ? undef : $INTEGRATED_conf_value;
                (doIt(-elem => "gerp"))   and $rH_data->{$key}->{GERP_RS_dbnsfp}            = ($GERP_RS_dbnsfp eq ".")            ? undef : $GERP_RS_dbnsfp;
                (doIt(-elem => "gerp"))   and $rH_data->{$key}->{GERP_rank_dbnsfp}          = ($GERP_rank_dbnsfp eq ".")          ? undef : $GERP_rank_dbnsfp;
                (doIt(-elem => "phylop")) and $rH_data->{$key}->{phyloP100way_vertebrate_score}     = ($phyloP100way_vertebrate_score eq ".")     ? undef : $phyloP100way_vertebrate_score;
                (doIt(-elem => "phylop")) and $rH_data->{$key}->{phyloP100way_vertebrate_rank}     = ($phyloP100way_vertebrate_rank eq ".")     ? undef : $phyloP100way_vertebrate_rank;
                (doIt(-elem => "phylop")) and $rH_data->{$key}->{phyloP20way_mammalian_score}     = ($phyloP20way_mammalian_score eq ".")     ? undef : $phyloP20way_mammalian_score;
                (doIt(-elem => "phylop")) and $rH_data->{$key}->{phyloP20way_mammalian_rank}     = ($phyloP20way_mammalian_rank eq ".")     ? undef : $phyloP20way_mammalian_rank;
		(doIt(-elem => "phcons")) and $rH_data->{$key}->{phastCons100way_vertebrate_score}  = ($phastCons100way_vertebrate_score eq ".")  ? undef : $phastCons100way_vertebrate_score;
		(doIt(-elem => "phcons")) and $rH_data->{$key}->{phastCons100way_vertebrate_rank}  = ($phastCons100way_vertebrate_rank eq "." ) ? undef : $phastCons100way_vertebrate_rank;
                (doIt(-elem => "phcons")) and $rH_data->{$key}->{phastCons20way_mammalian_score}  = ($phastCons20way_mammalian_score eq ".")  ? undef : $phastCons20way_mammalian_score;
                (doIt(-elem => "phcons")) and $rH_data->{$key}->{phastCons20way_mammalian_rank}  = ($phastCons20way_mammalian_rank eq "." ) ? undef : $phastCons20way_mammalian_rank;
                (doIt(-elem => "siphy"))  and $rH_data->{$key}->{SiPhy_29way_logOdds_score}= ($SiPhy_29way_logOdds_score eq ".")  ? undef : $SiPhy_29way_logOdds_score;
                (doIt(-elem => "siphy"))  and $rH_data->{$key}->{SiPhy_29way_logOdds_rank} = ($SiPhy_29way_logOdds_rank eq ".")   ? undef : $SiPhy_29way_logOdds_rank;
                (doIt(-elem => "interpro"))  and $rH_data->{$key}->{Interpro}              = ($Interpro eq ".")                   ? undef : $Interpro;
                (doIt(-elem => "gtex"))  and $rH_data->{$key}->{GTEx_V6_gene}              = ($GTEx_V6_gene eq ".")               ? undef : $GTEx_V6_gene;
                (doIt(-elem => "gtex"))  and $rH_data->{$key}->{GTEx_V6_tissue}            = ($GTEx_V6_tissue eq ".")             ? undef : $GTEx_V6_tissue;

            } else {

               # Dumper(%$rH_data->{$key});
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBNSFP);
    }
  
# DBNSFP end    
########################################################################################################################################################

    # if not excluded, run ANNOVAR -filter to output GERP++GT2 scores
    if (doIt(-elem => "gerp++gt2")) {
        my $annovarOutput = "$OUTFILE_GERP_GT2_WG.$USEDBUILD\_gerp++gt2_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GERP_GT2_WG', Number=>'A', Type=>'Float', Description=>'Whole-genome GERP++ scores greater than 2 (Annovar version 20140106)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype gerp++gt2
            $VERBOSE1 and printerr("\n**********\n* Resuming GERP++GT2 annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype gerp++gt2 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GERP_GT2_WG $ANNOINFILE $WG_GERPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GERP++GT2 scores with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GERP++GT2 predictions done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GERPGT2, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GERPGT2>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1         $2      $3    $4        $5     $6       $7         $8
            if ($line =~ /(gerp\+\+gt2)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GERP_GT2_WG} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GERPGT2);
    }

    # if not excluded, run ANNOVAR -filter to output ESP6500 scores
    if (doIt(-elem => "esp6500")) {
        my $annovarOutput = "$OUTFILE_ESP6500.$USEDBUILD\_esp6500si_all_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ESP6500', Number=>'A', Type=>'Float', Description=>'alternative allele frequency in all subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls (db updated 20141222)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype esp6500si_all
            $VERBOSE1 and printerr("\n**********\n* Resuming ESP6500 annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype esp6500si_all -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ESP6500 $ANNOINFILE $FREQUENCIESPATH/ESP/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ESP6500 frequencies with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of ESP6500 predictions done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ESP6500, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ESP6500>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1           $2     $3   $4         $5     $6       $7         $8
            if ($line =~ /(esp6500si_all)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ESP6500} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ESP6500);
    }

    # if not excluded, run ANNOVAR -filter to output ESP6500AA scores
    if (doIt(-elem => "esp6500aa")) {
        my $annovarOutput = "$OUTFILE_ESP6500AA.$USEDBUILD\_esp6500si_aa_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ESP6500AA', Number=>'A', Type=>'Float', Description=>'alternative allele frequency in African Americans in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls (db updated 20141222)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype esp6500si_aa
            $VERBOSE1 and printerr("\n**********\n* Resuming ESP6500AA annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype esp6500si_aa -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ESP6500AA $ANNOINFILE $FREQUENCIESPATH/ESP/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ESP6500AA frequencies with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of ESP6500AA predictions done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ESP6500AA, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ESP6500AA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1          $2     $3    $4        $5     $6       $7         $8
            if ($line =~ /(esp6500si_aa)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ESP6500AA} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ESP6500AA);
    }

    # if not excluded, run ANNOVAR -filter to output ESP6500EA scores
    if (doIt(-elem => "esp6500ea")) {
        my $annovarOutput = "$OUTFILE_ESP6500EA.$USEDBUILD\_esp6500si_ea_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ESP6500EA', Number=>'A', Type=>'Float', Description=>'alternative allele frequency in European Americans in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls (db updated 20141222)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype esp6500si_ea
            $VERBOSE1 and printerr("\n**********\n* Resuming ESP6500EA annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype esp6500si_ea -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ESP6500EA $ANNOINFILE $FREQUENCIESPATH/ESP/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ESP6500EA frequencies with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of ESP6500EA predictions done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ESP6500EA, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ESP6500EA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1          $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(esp6500si_ea)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ESP6500EA} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ESP6500EA);
    }

    # if not excluded, run ANNOVAR -filter to output GERP++ scores
    if (doIt(-elem => "gerp++elem")) {
        my $annovarOutput = "$OUTFILE_GERP_ELEM_WG.$USEDBUILD\_gerp++elem";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GERP_ELEM_WG', Number=>1, Type=>'Float', Description=>'Conserved genomic regions by GERP++ (Annovar version 20140106)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype gerp++elem
            $VERBOSE1 and printerr("\n**********\n* Resuming GERP++ELEM annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype gerp++elem -buildver $USEDBUILD -outfile $OUTFILE_GERP_ELEM_WG $ANNOINFILE $WG_GERPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GERP++ELEM conservation scores with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GERP++ELEM predictions done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GERPELEM, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GERPELEM>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1          $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(gerp\+\+elem)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                ($rH_data->{$key}->{GERP_ELEM_WG} = $2) =~ s/^Name\=//;;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GERPELEM);
    }

    # if not excluded, run ANNOVAR -filter to output REVEL annotations
    if (doIt(-elem => "revel")) {
        my $annovarOutput = "$OUTFILE_REVEL.$USEDBUILD\_revel_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'REVEL', Number=>1, Type=>'Float', Description=>'New ensemble method for predicting the pathogenicity of missense variants based on a combination of scores from 13 individual tools: MutPred, FATHMM v2.3, VEST 3.0, Polyphen-2, SIFT, PROVEAN, MutationAssessor, MutationTaster, LRT, GERP++, SiPhy, phyloP, and phastCons. version June-03-2016'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype revel
            $VERBOSE1 and printerr("\n**********\n* Resuming REVEL annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype revel -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_REVEL $ANNOINFILE $REVELPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of REVEL annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of REVEL annotations done...\n**********\n\n");
        }
        
        # open ANNOVAR output file for processing    
	open(ANNOVAR_REVEL, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_REVEL>) {
            chomp $line;    # $1     $2     $3     $4     $5     $6     $7     $8
            if ($line =~ /(revel)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{REVEL} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_REVEL);
    }

    # if not excluded, run ANNOVAR -filter to output Intervar annotations
    if (doIt(-elem => "intervar")) {
        my $annovarOutput = "$OUTFILE_INTERVAR.$USEDBUILD\_intervar_20170202_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'INTERVAR', Number=>'A', Type=>'String', Description=>'InterVar: clinical interpretation of missense variants (values: Uncertain significance, Benign, Likely benign, Likely pathogenic, Pathogenic). Annovar version 20170202'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype intervar
            $VERBOSE1 and printerr("\n**********\n* Resuming intervar annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype intervar_20170202 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_INTERVAR $ANNOINFILE $INTERVARPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of intervar annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of intervar annotations done...\n**********\n\n");
        }
        # open ANNOVAR output file for processing           
        open(ANNOVAR_INTERVAR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_INTERVAR>) {
            chomp $line;    # $1     $2     $3     $4     $5     $6     $7     $8
            if ($line =~ /(intervar_20170202)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{INTERVAR} = $2;
		$rH_data->{$key}->{INTERVAR} =~ s/ /_/g;
                $rH_data->{$key}->{INTERVAR} =~ s/__/_/g;
#                print $rH_data->{$key}->{INTERVAR} . "\n";
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_INTERVAR);
    }

    # if not excluded, run ANNOVAR -filter to output EBI3222PAK annotations
    if (doIt(-elem => "ebi3222pak")) {
        my $annovarOutput = "$OUTFILE_EBI3222PAK.$USEDBUILD\_ebi3222pak_20170902_dropped";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'EBI3222PAK', Number=>'A', Type=>'String', Description=>'Minor allele frequency of exome variants in 3222 Healthy Pakistani controls from EBI (exome captured with SureSelect_V5)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype ebi3222pak
            $VERBOSE1 and printerr("\n**********\n* Resuming ebi3222pak annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ebi3222pak_20170902 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EBI3222PAK $ANNOINFILE $FREQUENCIESPATH/EBI/ $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ebi3222pak annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of ebi3222pak annotations done...\n**********\n\n");
        }
        # open ANNOVAR output file for processing           
        open(ANNOVAR_EBI3222PAK, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_EBI3222PAK>) {
            chomp $line;    # $1     $2     $3     $4     $5     $6     $7     $8
            if ($line =~ /(ebi3222pak_20170902)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{EBI3222PAK} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_EBI3222PAK);
    }


    # if not excluded, run ANNOVAR -regionanno to output MIRNA regions
    if (doIt(-elem => "mirna")) {
        my $annovarOutput = "$OUTFILE_MIRNA.$USEDBUILD\_wgRna";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'MIRNA', Number=>1, Type=>'String', Description=>'snoRNA and miRNA annotations (from UCSC wgRna table, updated 2012Apr03)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype mirna
            $VERBOSE1 and printerr("\n**********\n* Resuming MIRNA annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirna -buildver $USEDBUILD -outfile $OUTFILE_MIRNA $ANNOINFILE $MIRNAPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of MIRNA annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of MIRNA annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_MIRNA, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_MIRNA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;  #  $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(mirna)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{MIRNA} = $2;
                $rH_data->{$key}->{MIRNA} =~ s/^Name\=//;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MIRNA);
    }

    # if not excluded, run ANNOVAR -regionanno to output MIRNATARGET regions
    if (doIt(-elem => "mirnatarget")) {
        my $annovarOutput = "$OUTFILE_MIRNATARGET.$USEDBUILD\_targetScanS";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'MIRNATARGET', Number=>1, Type=>'String', Description=>'TargetScan generated miRNA target site predictions (from UCSC TargetScanS table, updated 2012Apr03)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype mirnatarget
            $VERBOSE1 and printerr("\n**********\n* Resuming MIRNATARGET annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirnatarget -buildver $USEDBUILD -outfile $OUTFILE_MIRNATARGET $ANNOINFILE $MIRNAPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of MIRNATARGET annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of MIRNATARGET annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_MIRNATARGET, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_MIRNATARGET>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1         $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(mirnatarget)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{MIRNATARGET} = $2;
                $rH_data->{$key}->{MIRNATARGET} =~ s/^(Score\=\d+;)?Name\=//;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MIRNATARGET);
    }

    # if not excluded, run ANNOVAR -regionanno to output RNAseq regions from Cufflinks
    if (doIt(-elem => "rnaseq_cuff")) {
        my $annovarOutput = "$OUTFILE_RNASEQ_CUFF.$USEDBUILD\_bed";
        my $dbfile = "hg19_Cufflinks_allTissues.sorted.custom.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'RNASEQ_CUFF', Number=>1, Type=>'String', Description=>'Broad Cufflinks RNASeq alignment of BodyMap v2 (all tissues)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming RNASEQ_CUFF annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_CUFF $ANNOINFILE $BODYMAPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of RNASEQ_CUFF annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of RNASEQ_CUFF annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_RNASEQ_CUFF, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_RNASEQ_CUFF>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{RNASEQ_CUFF} = $2;
                $rH_data->{$key}->{RNASEQ_CUFF} =~ s/^Name\=//;
                ($rH_data->{$key}->{RNASEQ_CUFF} eq "NA") and $rH_data->{$key}->{RNASEQ_CUFF} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_RNASEQ_CUFF);
    }

    # if not excluded, run ANNOVAR -regionanno to output RNAseq regions from Scripture
    if (doIt(-elem => "rnaseq_scri")) {
        my $annovarOutput = "$OUTFILE_RNASEQ_SCRI.$USEDBUILD\_bed";
        my $dbfile = "hg19_Scripture_allTissues.sorted.custom.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'RNASEQ_SCRI', Number=>1, Type=>'String', Description=>'Broad Scripture RNASeq alignment of BodyMap v2 (all tissues)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming RNASEQ_SCRI annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_SCRI $ANNOINFILE $BODYMAPPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of RNASEQ_SCRI annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of RNASEQ_SCRI annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_RNASEQ_SCRI, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_RNASEQ_SCRI>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3   $4         $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{RNASEQ_SCRI} = $2;
                $rH_data->{$key}->{RNASEQ_SCRI} =~ s/^Name\=//;
                ($rH_data->{$key}->{RNASEQ_SCRI} eq "NA") and $rH_data->{$key}->{RNASEQ_SCRI} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_RNASEQ_SCRI);
    }

    # if not excluded, run ANNOVAR -regionanno to output CLINVAR Description
    if (doIt(-elem => "clinvar")) {
        my $annovarOutput = "$OUTFILE_CLINVAR.$USEDBUILD\_clinvar_20170130_dropped";
        my $dbfile = "hg19_clinvar_20170130.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'CLINVAR', Number=>'A', Type=>'String', Description=>'Clinically significant variants from NCBI\'s ClinVar resource (ANNOVAR - clinvar_20170130) with infos CLINSIG CLNDBN  CLNACC  CLNDSDB CLNDSDBID'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype clinvar_20170130
            $VERBOSE1 and printerr("\n**********\n* Resuming CLINVAR annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype clinvar_20170130 -otherinfo -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_CLINVAR $ANNOINFILE $CLINVARPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of CLINVAR annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of CLINVAR annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_CLINVAR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_CLINVAR>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1                $2    $3    $4        $5     $6       $7           $8
            if ($line =~ /(clinvar_20170130)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{CLINVAR} = $2;
                $rH_data->{$key}->{CLINVAR} =~ s/\=/\:/g;
                $rH_data->{$key}->{CLINVAR} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{CLINVAR} =~ s/ /_/g;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_CLINVAR);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax ChipSeq Description
    if (doIt(-elem => "gtx_chipseq")) {
        my $annovarOutput = "$OUTFILE_GTX_CHIPSEQ.$USEDBUILD\_bed";
        my $dbfile = "chip_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_CHIPSEQ', Number=>1, Type=>'String', Description=>'Genome Trax Predicted ChIP-Seq TFBS description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_CHIPSEQ annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CHIPSEQ $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_CHIPSEQ annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_CHIPSEQ annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_CHIPSEQ, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_CHIPSEQ>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_CHIPSEQ} = $2;
                $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_CHIPSEQ} eq "NA") and $rH_data->{$key}->{GTX_CHIPSEQ} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_CHIPSEQ);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax CPG Description
    if (doIt(-elem => "gtx_cpg")) {
        my $annovarOutput = "$OUTFILE_GTX_CPG.$USEDBUILD\_bed";
        my $dbfile = "cpg_islands_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_CPG', Number=>1, Type=>'String', Description=>'Genome Trax CpG Islands description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_CPG annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CPG $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_CPG annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_CPG annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_CPG, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_CPG>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_CPG} = $2;
                $rH_data->{$key}->{GTX_CPG} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_CPG} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_CPG} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_CPG} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_CPG} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_CPG} eq "NA") and $rH_data->{$key}->{GTX_CPG} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_CPG);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax Disease Description
    if (doIt(-elem => "gtx_disease")) {
        my $annovarOutput = "$OUTFILE_GTX_DISEASE.$USEDBUILD\_bed";
        my $dbfile = "disease_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_DISEASE', Number=>1, Type=>'String', Description=>'Genome Trax Disease associations description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_DISEASE annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DISEASE $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_DISEASE annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_DISEASE annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_DISEASE, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_DISEASE>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_DISEASE} = $2;
                $rH_data->{$key}->{GTX_DISEASE} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_DISEASE} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_DISEASE} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_DISEASE} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_DISEASE} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_DISEASE} eq "NA") and $rH_data->{$key}->{GTX_DISEASE} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_DISEASE);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax Dnase binding sites
    if (doIt(-elem => "gtx_dnase")) {
        my $annovarOutput = "$OUTFILE_GTX_DNASE.$USEDBUILD\_bed";
        my $dbfile = "dnase_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_DNASE', Number=>1, Type=>'String', Description=>'Genome Trax predicted binding sites at the DNAse hypersensitive regions (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_DNASE annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DNASE $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_DNASE annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_DNASE annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_DNASE, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_DNASE>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_DNASE} = $2;
                $rH_data->{$key}->{GTX_DNASE} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_DNASE} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_DNASE} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_DNASE} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_DNASE} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_DNASE} eq "NA") and $rH_data->{$key}->{GTX_DNASE} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_DNASE);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax Drug Description
    if (doIt(-elem => "gtx_drug")) {
        my $annovarOutput = "$OUTFILE_GTX_DRUG.$USEDBUILD\_bed";
        my $dbfile = "drug_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_DRUG', Number=>1, Type=>'String', Description=>'Genome Trax Drug Targets description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_DRUG annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DRUG $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_DRUG annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_DRUG annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_DRUG, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_DRUG>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_DRUG} = $2;
                $rH_data->{$key}->{GTX_DRUG} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_DRUG} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_DRUG} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_DRUG} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_DRUG} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_DRUG} eq "NA") and $rH_data->{$key}->{GTX_DRUG} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_DRUG);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax GWAS Description
    if (doIt(-elem => "gtx_gwas")) {
        my $annovarOutput = "$OUTFILE_GTX_GWAS.$USEDBUILD\_generic_dropped";
        my $dbfile = "gwas_hg19.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_GWAS', Number=>'1', Type=>'String', Description=>'Genome Trax GWAS Catalogue description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_GWAS annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GTX_GWAS $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_GWAS annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_GWAS annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_GWAS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_GWAS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;  #  $1       $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(generic)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_GWAS} = $2;
                $rH_data->{$key}->{GTX_GWAS} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_GWAS} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_GWAS} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_GWAS} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_GWAS} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_GWAS} eq "NA") and $rH_data->{$key}->{GTX_GWAS} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_GWAS);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_hgmd_disgene")) {
        my $annovarOutput = "$OUTFILE_GTX_HGMD_DISGENE.$USEDBUILD\_bed";
        my $dbfile = "hgmd_disease_genes_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_HGMD_DISGENE', Number=>1, Type=>'String', Description=>'Genome Trax HGMD disease genes description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_HGMD_DISGENE annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD_DISGENE $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_HGMD_DISGENE annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_HGMD_DISGENE annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_HGMD_DISGENE, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_HGMD_DISGENE>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_HGMD_DISGENE} = $2;
                $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_HGMD_DISGENE} eq "NA") and $rH_data->{$key}->{GTX_HGMD_DISGENE} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_HGMD_DISGENE);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMDIMPUTED Description
    if (doIt(-elem => "gtx_hgmdimputed")) {
        my $annovarOutput = "$OUTFILE_GTX_HGMDIMPUTED.$USEDBUILD\_generic_dropped";
        my $dbfile = "hgmdimputed_hg19.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_HGMDIMPUTED', Number=>'1', Type=>'String', Description=>'Genome Trax HGMDIMPUTED calculate all alternative possible codon changes that result in the same amino acid change reported in the actual HGMD record. To each of these imputed mutations, the information from the corresponding mutation from HGMD is transferred (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_HGMDIMPUTED annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GTX_HGMDIMPUTED $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_HGMDIMPUTED annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_HGMDIMPUTED annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_HGMDIMPUTED, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_HGMDIMPUTED>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;  #  $1       $2     $3     $4       $5     $6       $7           $8
            if ($line =~ /(generic)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_HGMDIMPUTED} = $2;
                $rH_data->{$key}->{GTX_HGMDIMPUTED} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_HGMDIMPUTED} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_HGMDIMPUTED} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_HGMDIMPUTED} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_HGMDIMPUTED} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_HGMDIMPUTED} eq "NA") and $rH_data->{$key}->{GTX_HGMDIMPUTED} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_HGMDIMPUTED);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Description
    if (doIt(-elem => "gtx_hgmd")) {
        my $annovarOutput = "$OUTFILE_GTX_HGMD.$USEDBUILD\_generic_dropped";
        my $dbfile = "hgmd_hg19.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_HGMD', Number=>'A', Type=>'String', Description=>'Genome Trax HGMD inherited (germ-line) disease mutations description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_HGMD annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GTX_HGMD $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_HGMD annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_HGMD annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_HGMD, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_HGMD>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;  #  $1       $2     $3     $4       $5     $6       $7           $8
            if ($line =~ /(generic)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_HGMD} = $2;
                $rH_data->{$key}->{GTX_HGMD} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_HGMD} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_HGMD} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_HGMD} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_HGMD} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_HGMD} eq "NA") and $rH_data->{$key}->{GTX_HGMD} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_HGMD);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax Microsatellite Description
    if (doIt(-elem => "gtx_microsat")) {
        my $annovarOutput = "$OUTFILE_GTX_MICROSAT.$USEDBUILD\_bed";
        my $dbfile = "microsatellites_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_MICROSAT', Number=>1, Type=>'String', Description=>'Genome Trax Microsatellite repeats description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_MICROSAT annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MICROSAT $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_MICROSAT annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_MICROSAT annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_MICROSAT, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_MICROSAT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_MICROSAT} = ($2) ? 1 : undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_MICROSAT);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax miRNA Description
    if (doIt(-elem => "gtx_mirna")) {
        my $annovarOutput = "$OUTFILE_GTX_MIRNA.$USEDBUILD\_bed";
        my $dbfile = "mirna_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_MIRNA', Number=>1, Type=>'String', Description=>'Genome Trax microRNA sequences description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_MIRNA annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MIRNA $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_MIRNA annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_MIRNA annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_MIRNA, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_MIRNA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_MIRNA} = $2;
                $rH_data->{$key}->{GTX_MIRNA} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_MIRNA} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_MIRNA} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_MIRNA} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_MIRNA} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_MIRNA} eq "NA") and $rH_data->{$key}->{GTX_MIRNA} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_MIRNA);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax OMIM Description
    if (doIt(-elem => "gtx_omim")) {
        my $annovarOutput = "$OUTFILE_GTX_OMIM.$USEDBUILD\_bed";
        my $dbfile = "omim_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_OMIM', Number=>1, Type=>'String', Description=>'Genome Trax OMIM disorders (BIOBASE - 15th September 2015)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_OMIM annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_OMIM $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_OMIM annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_OMIM annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_OMIM, "$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_OMIM>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_OMIM} = $2;
                $rH_data->{$key}->{GTX_OMIM} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_OMIM} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_OMIM} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_OMIM} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_OMIM} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_OMIM} eq "NA") and $rH_data->{$key}->{GTX_OMIM} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_OMIM);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax ORPHA Description
    
    if (doIt(-elem => "gtx_orpha")) {
        my $annovarOutput = "$OUTFILE_GTX_ORPHA.$USEDBUILD\_bed";
        my $dbfile = "orpha_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_ORPHA', Number=>1, Type=>'String', Description=>'Genome Trax ORPHA information on rare diseases and orphan drugs (BIOBASE - 15th September 2015)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_ORPHA annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_ORPHA $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_ORPHA annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_ORPHA annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_ORPHA, "$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_ORPHA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_ORPHA} = $2;
                $rH_data->{$key}->{GTX_ORPHA} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_ORPHA} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_ORPHA} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_ORPHA} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_ORPHA} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_ORPHA} eq "NA") and $rH_data->{$key}->{GTX_ORPHA} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_ORPHA);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax Pathway Description
    if (doIt(-elem => "gtx_path")) {
        my $annovarOutput = "$OUTFILE_GTX_PATH.$USEDBUILD\_bed";
        my $dbfile = "pathway_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_PATH', Number=>1, Type=>'String', Description=>'Genome Trax Pathways membership description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_PATH annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PATH $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_PATH annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_PATH annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_PATH, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_PATH>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_PATH} = $2;
                $rH_data->{$key}->{GTX_PATH} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_PATH} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_PATH} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_PATH} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_PATH} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_PATH} eq "NA") and $rH_data->{$key}->{GTX_PATH} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_PATH);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax PGMD Description
    if (doIt(-elem => "gtx_pgmd")) {
        my $annovarOutput = "$OUTFILE_GTX_PGMD.$USEDBUILD\_bed";
        my $dbfile = "pgmd_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_PGMD', Number=>1, Type=>'String', Description=>'Genome Trax PGMD indicate variants that have been shown to exhibit a pharmacogenomic effect on patients (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_PGMD annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PGMD $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_PGMD annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) { 
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_PGMD annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_PGMD, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_PGMD>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_PGMD} = $2;
                $rH_data->{$key}->{GTX_PGMD} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_PGMD} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_PGMD} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_PGMD} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_PGMD} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_PGMD} eq "NA") and $rH_data->{$key}->{GTX_PGMD} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_PGMD);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax PTMS Description
    if (doIt(-elem => "gtx_ptms")) {
        my $annovarOutput = "$OUTFILE_GTX_PTMS.$USEDBUILD\_bed";
        my $dbfile = "ptms_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_PTMS', Number=>1, Type=>'String', Description=>'Genome Trax Post translational modifications description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_PTMS annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PTMS $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_PTMS annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_PTMS annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_PTMS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_PTMS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_PTMS} = $2;
                $rH_data->{$key}->{GTX_PTMS} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_PTMS} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_PTMS} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_PTMS} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_PTMS} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_PTMS} eq "NA") and $rH_data->{$key}->{GTX_PTMS} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_PTMS);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax Transfac Description
    if (doIt(-elem => "gtx_transfac")) {
        my $annovarOutput = "$OUTFILE_GTX_TRANSFAC_TFBS.$USEDBUILD\_bed";
        my $dbfile = "transfac_sites_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_TRANSFAC_TFBS', Number=>1, Type=>'String', Description=>'Genome Trax TRANSFAC experimentally verified TFBS description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_TRANSFAC_TFBS annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TRANSFAC_TFBS $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_TRANSFAC_TFBS annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_TRANSFAC_TFBS annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_TRANSFAC_TFBS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_TRANSFAC_TFBS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_TRANSFAC_TFBS} = $2;
                $rH_data->{$key}->{GTX_TRANSFAC_TFBS} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_TRANSFAC_TFBS} =~ s/\s+;/\|/g;
                $rH_data->{$key}->{GTX_TRANSFAC_TFBS} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_TRANSFAC_TFBS} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_TRANSFAC_TFBS} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_TRANSFAC_TFBS} eq "NA") and $rH_data->{$key}->{GTX_TRANSFAC_TFBS} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_TRANSFAC_TFBS);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax TSS Description
    if (doIt(-elem => "gtx_tss")) {
        my $annovarOutput = "$OUTFILE_GTX_TSS.$USEDBUILD\_bed";
        my $dbfile = "tss_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_TSS', Number=>1, Type=>'String', Description=>'Genome Trax TSSs (transcription start sites) description (BIOBASE - GenomeTrax.2017.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_TSS annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TSS $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_TSS annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_TSS annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_TSS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_TSS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_TSS} = $2;
                $rH_data->{$key}->{GTX_TSS} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_TSS} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_TSS} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_TSS} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_TSS} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_TSS} eq "NA") and $rH_data->{$key}->{GTX_TSS} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_TSS);
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax VISTA Description
    if (doIt(-elem => "gtx_vista")) {
        my $annovarOutput = "$OUTFILE_GTX_VISTA.$USEDBUILD\_bed";
        my $dbfile = "vista_hg19.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_VISTA', Number=>1, Type=>'String', Description=>'Genome Trax VISTA track helps predict the damaging effect of a variation occurring in an enhancer site that in turn could affect gene expression (BIOBASE - GenomeTrax.2015.2)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming GTX_VISTA annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_VISTA $ANNOINFILE $GENOMETRAXPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of GTX_VISTA annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of GTX_VISTA annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_GTX_VISTA, "$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_GTX_VISTA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GTX_VISTA} = $2;
                $rH_data->{$key}->{GTX_VISTA} =~ s/^Name\=//;
                $rH_data->{$key}->{GTX_VISTA} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{GTX_VISTA} =~ s/, /,/g;
                $rH_data->{$key}->{GTX_VISTA} =~ s/ /_/g;
                $rH_data->{$key}->{GTX_VISTA} =~ s/__/_/g;
                ($rH_data->{$key}->{GTX_VISTA} eq "NA") and $rH_data->{$key}->{GTX_VISTA} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GTX_VISTA);
    }


    # if not excluded, run ANNOVAR -regionanno to output Encode Chromatin States Description
    if (doIt(-elem => "enc_hmm")) {
        my $annovarOutput = "$OUTFILE_ENC_HMM.$USEDBUILD\_bed";
        my $dbfile = "hg19_wgEncodeBroadHmm_NoCancerCombined_severeStates_sorted.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ENC_HMM', Number=>1, Type=>'String', Description=>'Encode top 5 Chromatin States from non cancer cell lines : concatenation of state & cell line (ENCODE - wgEncodeBroadHmm[Gm12878|H1hesc|Huvec|Hsmm|Huvec|Nhek|Nhlf]HMM tracks)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming ENC_HMM annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_ENC_HMM $ANNOINFILE $ENCODEPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ENC_HMM annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of ENC_HMM annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ENC_HMM, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ENC_HMM>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $annoString = $2;
                $annoString =~ s/^Name\=//;
                $annoString =~ s/\s*;/\|/g;
                $annoString =~ s/, /,/g;
                $annoString =~ s/ /_/g;
                $annoString =~ s/__/_/g;
                if ($annoString eq "NA") {
                    $rH_data->{$key}->{ENC_HMM} = undef;
                } else {
                    my @annots = split(/,/, $annoString);
                    my %annoHash = ();
                    foreach my $annot (@annots) {
                        $annot =~ /(.+)\((.+)\)/;
                        (defined $annoHash{$1}) and $annoHash{$1} .= "," . $2;
                        (!defined $annoHash{$1}) and $annoHash{$1} = $2;
                    }
                    $rH_data->{$key}->{ENC_HMM} = join ",", map { $_ . "(" . $annoHash{$_} . ")" } sort keys %annoHash;
                }
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ENC_HMM);
    }

    # if not excluded, run ANNOVAR -regionanno to output Encode Transcription Factor Binding Sites Description
    if (doIt(-elem => "enc_tfbs")) {
        my $annovarOutput = "$OUTFILE_ENC_TFBS.$USEDBUILD\_bed";
        my $dbfile = "hg19_wgEncodeRegTfbsClusteredWithCellsV3_formated.bed";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ENC_TFBS', Number=>1, Type=>'String', Description=>'Transcription factor binding sites using ChIP-seq from the ENCODE project : concatenation of item name, score & cell line(s) (ENCODE - wgEncodeRegTfbsClusteredWithCellsV3 track)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming ENC_TFBS annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_ENC_TFBS $ANNOINFILE $ENCODEPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of ENC_TFBS annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of ENC_TFBS annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_ENC_TFBS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_ENC_TFBS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ENC_TFBS} = $2;
                $rH_data->{$key}->{ENC_TFBS} =~ s/^Name\=//;
                $rH_data->{$key}->{ENC_TFBS} =~ s/\s*;/\|/g;
                $rH_data->{$key}->{ENC_TFBS} =~ s/, /,/g;
                $rH_data->{$key}->{ENC_TFBS} =~ s/ /_/g;
                $rH_data->{$key}->{ENC_TFBS} =~ s/__/_/g;
                ($rH_data->{$key}->{ENC_TFBS} eq "NA") and $rH_data->{$key}->{ENC_TFBS} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ENC_TFBS);
    }

    # if not excluded, run ANNOVAR -regionanno to output Encode Transcription Factor Binding Sites Description
    if (doIt(-elem => "simple_repeat")) {
        my $annovarOutput = "$OUTFILE_SIMPLE_REPEAT.$USEDBUILD\_simpleRepeat";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'SIMPLE_REPEAT', Number=>1, Type=>'String', Description=>'Simple Tandem Repeats located by Tandem Repeats Finder (TRF) : concatenation of percent match, percent indel & repeat sequence (UCSC - simpleRepeats track)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype simpleRepeat
            $VERBOSE1 and printerr("\n**********\n* Resuming SIMPLE_REPEAT annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype simpleRepeat -colsWanted 8,6,16 -buildver $USEDBUILD -outfile $OUTFILE_SIMPLE_REPEAT $ANNOINFILE $REPEATPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of SIMPLE_REPEAT annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of SIMPLE_REPEAT annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_SIMPLE_REPEAT, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_SIMPLE_REPEAT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; #      $1         $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(simpleRepeat)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{SIMPLE_REPEAT} = $2;
                $rH_data->{$key}->{SIMPLE_REPEAT} =~ s/^Name\=(\d+)/$1\%/;
                $rH_data->{$key}->{SIMPLE_REPEAT} =~ s/,(\d+)\:(\d+\.?\d*)\:(\w+)/,$1\%\:$2\:$3/g;
                ($rH_data->{$key}->{SIMPLE_REPEAT} eq "NA") and $rH_data->{$key}->{SIMPLE_REPEAT} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_SIMPLE_REPEAT);
    }

    # if not excluded, run ANNOVAR -regionanno to output Encode Transcription Factor Binding Sites Description
    if (doIt(-elem => "segment_dups")) {
        my $annovarOutput = "$OUTFILE_SEGMENTAL_DUPS.$USEDBUILD\_genomicSuperDups";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'SEGMENT_DUPS', Number=>1, Type=>'String', Description=>'Regions detected as putative genomic duplications : concatenation of fraction match & genomic location of the duplication (UCSC - genomicSuperDups track)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype genomicSuperDups
            $VERBOSE1 and printerr("\n**********\n* Resuming SEGMENT_DUPS annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype genomicSuperDups -colsWanted 7,8,9 -scorecolumn 26 -buildver $USEDBUILD -outfile $OUTFILE_SEGMENTAL_DUPS $ANNOINFILE $REPEATPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of SEGMENT_DUPS annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of SEGMENT_DUPS annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_SEGMENT_DUPS, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_SEGMENT_DUPS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #     $1           $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(genomicSuperDups)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{SEGMENT_DUPS} = $2;
                $rH_data->{$key}->{SEGMENT_DUPS} =~ s/^Score\=//;
                $rH_data->{$key}->{SEGMENT_DUPS} =~ s/Name\=//;
                $rH_data->{$key}->{SEGMENT_DUPS} =~ s/\s*;/\:/g;
                $rH_data->{$key}->{SEGMENT_DUPS} =~ s/, /,/g;
                $rH_data->{$key}->{SEGMENT_DUPS} =~ s/ /_/g;
                $rH_data->{$key}->{SEGMENT_DUPS} =~ s/__/_/g;
                $rH_data->{$key}->{SEGMENT_DUPS} =~ s/(chr\d+):(\d+):(\d+)/$1\_$2\_$3/g;
                ($rH_data->{$key}->{SEGMENT_DUPS} eq "NA") and $rH_data->{$key}->{SEGMENT_DUPS} = undef;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_SEGMENT_DUPS);
    }

    # if not excluded, run ANNOVAR -regionanno to output Encode Transcription Factor Binding Sites Description
    if (doIt(-elem => "rmsk")) {
        my $annovarOutput = "$OUTFILE_REPEAT_MASKER.$USEDBUILD\_rmsk_custom";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'REPEAT_MASKER', Number=>1, Type=>'String', Description=>'Regions detected by RepeatMasker as interspersed repeats and low complexity DNA sequences : concatenation of the repeat Name, Class & Family, followed by its genomic location (UCSC - rmsk track)'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype rmsk_custom
            $VERBOSE1 and printerr("\n**********\n* Resuming REPEAT_MASKER annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype rmsk_custom -buildver $USEDBUILD -outfile $OUTFILE_REPEAT_MASKER $ANNOINFILE $REPEATPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of REPEAT_MASKER annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of REPEAT_MASKER annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_REPEAT_MASKER, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_REPEAT_MASKER>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #   $1         $2    $3    $4        $5     $6       $7           $8
            if ($line =~ /(rmsk_custom)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{REPEAT_MASKER} = $2;
                $rH_data->{$key}->{REPEAT_MASKER} =~ s/Name\=//;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_REPEAT_MASKER);
    }

    # if not excluded, run ANNOVAR -regionanno to output RVIS ExAC score Description
    if (doIt(-elem => "rvis_exac")) {
        my $annovarOutput = "$OUTFILE_RVIS_EXAC.$USEDBUILD\_bed";
        my $dbfile = "GenicIntolerance_v3_12Mar16.annovar_ready.filtered.split.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'RVIS_raw', Number=>1, Type=>'Float', Description=>'RVIS Raw score - version 3 12Mar16'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'RVIS_percent', Number=>1, Type=>'Float', Description=>'RVIS Percentile rank - version 3 12Mar16'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_RVIS_percent', Number=>1, Type=>'Float', Description=>'RVIS from ExAC score - version 3 12Mar16'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ExAC_LoF_FDR', Number=>1, Type=>'Float', Description=>'ExAC LoF FDR adjusted p-value reflects the significance of the departure from the expected rate of LoF variants - version 3 12Mar16'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -regionanno -dbtype bed
            $VERBOSE1 and printerr("\n**********\n* Resuming RVIS_ExAC annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RVIS_EXAC $ANNOINFILE $RVIS_EXACPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of RVIS_ExAC annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of RVIS_ExAC annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_RVIS_EXAC, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_RVIS_EXAC>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; # $1     $2     $3    $4        $5     $6       $7           $8
            if ($line =~ /(bed)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $rvis_annot = $2;
                $rvis_annot =~ s/^Name\=//;
                
                my ($RVIS_raw, $RVIS_percent, $ExAC_RVIS_percent, $ExAC_LoF_FDR) = split(":",$rvis_annot);
                $rH_data->{$key}->{RVIS_raw} = $RVIS_raw;
                ($rH_data->{$key}->{RVIS_raw} eq "NA") and $rH_data->{$key}->{RVIS_raw} = undef;
                $rH_data->{$key}->{RVIS_percent} = $RVIS_percent;
                ($rH_data->{$key}->{RVIS_percent} eq "NA") and $rH_data->{$key}->{RVIS_percent} = undef;
                $rH_data->{$key}->{ExAC_RVIS_percent} = $ExAC_RVIS_percent;
                ($rH_data->{$key}->{ExAC_RVIS_percent} eq "NA") and $rH_data->{$key}->{ExAC_RVIS_percent} = undef;
                $rH_data->{$key}->{ExAC_LoF_FDR} = $ExAC_LoF_FDR;
                ($rH_data->{$key}->{ExAC_LoF_FDR} eq "N/A") and $rH_data->{$key}->{ExAC_LoF_FDR} = undef;

            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_RVIS_EXAC);
    }

    # if not excluded, run ANNOVAR -filter to output MPC ExAC score Description
    if (doIt(-elem => "mpc_exac")) {
        my $annovarOutput = "$OUTFILE_MPC_EXAC.$USEDBUILD\_generic_dropped";
        my $dbfile = "exac_regional_missense_constraint_v2.mpc_values.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'MPC_EXAC', Number=>1, Type=>'Float', Description=>'Exac regional missense constraint transformed fitted MPC score - version 2 14Jul17'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming MPC_ExAC annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -otherinfo -buildver $USEDBUILD -outfile $OUTFILE_MPC_EXAC $ANNOINFILE $MPC_EXACPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of MPC_ExAC annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of MPC_ExAC annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_MPC_EXAC, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_MPC_EXAC>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; #   $1       $2     $3     $4         $5     $6       $7           $8
            if ($line =~ /(generic)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $MPC_EXAC= $2;
                
                $rH_data->{$key}->{MPC_EXAC} = $MPC_EXAC;
                ($rH_data->{$key}->{MPC_EXAC} eq "NA") and $rH_data->{$key}->{MPC_EXAC} = undef;

            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MPC_EXAC);
    }

    # if not excluded, run ANNOVAR -filter to output uORF scores
    if (doIt(-elem => "uorf")) {
        my $annovarOutput = "$OUTFILE_UORF.$USEDBUILD\_generic_dropped";
        my $dbfile = "uORF.annovar.generic_db.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'uORF_type', Number=>'A', Type=>'String', Description=>'5\' UTR variant impact from Whiffin et al. 2019: effect type (stop vs uAUG-creating)' });
        $rO_vcf->add_header_line({key=>'INFO', ID=>'uORF_geneClass', Number=>'A', Type=>'String', Description=>'5\' UTR variant impact from Whiffin et al. 2019: gene class'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'uORF_distance', Number=>'A', Type=>'Float', Description=>'5\' UTR variant impact from Whiffin et al. 2019: distance to Stop/CDS'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'uORF_effect', Number=>'A', Type=>'String', Description=>'5\' UTR variant impact from Whiffin et al. 2019: effect'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'uORF_kozakStrength', Number=>'A', Type=>'String', Description=>'5\' UTR variant impact from Whiffin et al. 2019: Kozak Strength'});
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming uORF annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -otherinfo -buildver $USEDBUILD -outfile $OUTFILE_UORF $ANNOINFILE $UORFPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of uORF annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of uORF annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_UORF, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_UORF>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line; #   $1       $2     $3     $4         $5     $6       $7           $8
            if ($line =~ /(generic)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $UORF= $2;
                my ($type, $class, $distance, $effect, $strength) = split(':',$UORF);

                
                $rH_data->{$key}->{uORF_type} = $type;
                $rH_data->{$key}->{uORF_geneClass} = $class;
                $rH_data->{$key}->{uORF_distance} = $distance;
                $rH_data->{$key}->{uORF_effect} = $effect;
                $rH_data->{$key}->{uORF_kozakStrength} = $strength;
                ($rH_data->{$key}->{uORF_distance} eq "NA") and $rH_data->{$key}->{uORF_distance} = undef;

            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_UORF);
    }

    # if not excluded, run ANNOVAR -filter to output UV_signature mutations
    if (doIt(-elem => "uv_signature")) {
        my $annovarOutput = "$OUTFILE_UV_SIGNATURE.$USEDBUILD\_generic_dropped";
        my $dbfile = "v37.UV_signature.annovar.generic_db.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'UV_signature', Number=>'A', Type=>'String', Description=>'Dipyrimidine transversion point mutations' });
        if ($RESUMEFOLDER && -e $annovarOutput) {
            # resuming ANNOVAR -filter -dbtype generic
            $VERBOSE1 and printerr("\n**********\n* Resuming UV_signature annotation using already existing file $annovarOutput\n**********\n\n");
        } else {
            my $systemCommand = "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -otherinfo -buildver $USEDBUILD -outfile $OUTFILE_UV_SIGNATURE $ANNOINFILE $UVSIGNATUREPATH $VERBOSE1";
            $VERBOSE1 and printerr("\n**********\n* Starting fetching of UV_signature annotations with system command:\n* <$systemCommand>\n**********\n\n");
            my $error = system ($systemCommand);
            if ($error) {
                printerr("\n**********\n* Error (code: $error) while running system command: <$systemCommand>\n$!\n**********\n\n");
            } elsif ($VERBOSE2) {
                my $exitSignal = ($? >> 8);
                my $message = sprintf "Command exited with value %d", $exitSignal;
                printerr("\n**********\n* $message\n**********\n\n");
            }
            $VERBOSE1 and printerr("\n**********\n* Fetching of UV_signature annotations done...\n**********\n\n");
        }

        # open ANNOVAR output file for processing
        open(ANNOVAR_UV_SIGNATURE, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        $VERBOSE1 and printerr("\n**********\n* Parsing $annovarOutput to populate variant hash...\n**********\n\n");
        while (my $line = <ANNOVAR_UV_SIGNATURE>) {
            chomp $line; #   $1       $2     $3     $4         $5     $6       $7           $8
            if ($line =~ /(generic)\t(.+?)\t(.+?([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $UV_SIGNATURE= $2;
                
                $rH_data->{$key}->{UV_signature} = $UV_SIGNATURE;

            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_UV_SIGNATURE);
    }

}

sub printVcf {
    my %arg = @_;
    my $rO_vcf      = $arg{-vcf};
    my $rH_data     = $arg{-data};
    my $rH_vcfIndex = $arg{-index};

    # output results to a results file
    my $resultFile = $OUTFILE;
    (not $resultFile) and $resultFile = ($VARIANT_CALLER) ? $OUTPATH . $ID . "_ANNOVAR_output." . $VARIANT_CALLER . "." . $BUILDVER . ".RESULTS.vcf" : $OUTPATH . $ID . "_ANNOVAR_output." . $BUILDVER . ".RESULTS.vcf";
    open(RESULTS, ">$resultFile") || die ("Could not open RESULTS file $resultFile :\n$!\n");

    print RESULTS $rO_vcf->format_header();

    my $i = 0;
    while (my $rA_x = $rO_vcf->next_data_array()) {
        $i++;
        my ($VC, $VT, $VFT, $Gene, $DA, $oneKG, $CG46, $NCI60, $COSMIC, $ICGC21_ID, $ICGC21_Occurrence, $GoNL, $CLINVAR) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($WGSCADD, $WGSCADD_Phred, $WGSDANN, $dbscSNV_ADA_SCORE, $dbscSNV_RF_SCORE, $HRCR1_AF, $HRCR1_AC, $HRCR1_AN, $HRCR1_nonKG_AF, $HRCR1_nonKG_AC, $HRCR1_nonKG_AN , $WGSFATHMM_coding, $WGSFATHMM_noncoding, $GWAVA_region_score, $GWAVA_tss_score,$ExAC_ALL_nonpsych, $ExAC_ALL_nontcga, $KAVIAR_AF, $KAVIAR_AC, $KAVIAR_AN, $WGSEIGEN) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($ExAC_Freq, $ExAC_AFR, $ExAC_AMR, $ExAC_EAS, $ExAC_FIN, $ExAC_NFE, $ExAC_OTH, $ExAC_SAS) = (undef, undef, undef, undef, undef, undef, undef, undef);
        my ($GME_AF, $GME_NWA, $GME_NEA, $GME_AP,  $GME_Israel, $GME_SD, $GME_TP, $GME_CA) = (undef, undef, undef, undef, undef, undef, undef, undef);
        my ($AFR, $AMR, $EUR, $EAS, $SAS, $RNASEQ_CUFF, $RNASEQ_SCRI, $SIFT_pred, $SIFT_rank, $SIFT_score) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($dbSNP_132, $dbSNP_138, $dbSNP_138Common, $dbSNP_138CommonAlCount, $dbSNP_138CommonAlFreq, $dbSNP_138Flagged) = (undef, undef, undef, undef, undef, undef);
        my ($dbSNP_138FlaggedAlCount, $dbSNP_138FlaggedAlFreq, $ESP6500, $ESP6500AA, $ESP6500EA, $MIRNA, $MIRNATARGET) = (undef, undef, undef, undef, undef, undef, undef);
        my ($dbSNP_142, $dbSNP_142Common, $dbSNP_142CommonAlCount, $dbSNP_142CommonAlFreq, $dbSNP_142Flagged,$dbSNP_142FlaggedAlCount, $dbSNP_142FlaggedAlFreq) = (undef, undef, undef, undef, undef, undef, undef);
        my ($gnomAD_exome_ALL, $gnomAD_exome_AFR, $gnomAD_exome_AMR, $gnomAD_exome_ASJ, $gnomAD_exome_EAS, $gnomAD_exome_FIN, $gnomAD_exome_NFE, $gnomAD_exome_OTH, $gnomAD_exome_SAS, $gnomAD_genome_ALL, $gnomAD_genome_AFR, $gnomAD_genome_AMR, $gnomAD_genome_ASJ, $gnomAD_genome_EAS, $gnomAD_genome_FIN, $gnomAD_genome_NFE, $gnomAD_genome_OTH) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($Polyphen2_HDIV_pred, $Polyphen2_HDIV_score, $Polyphen2_HDIV_rank, $Polyphen2_HVAR_pred, $Polyphen2_HVAR_score, $Polyphen2_HVAR_rank, $LRT_pred, $LRT_score, $LRT_rank) = (undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($MutationTaster_pred, $MutationTaster_score, $MutationTaster_rank, $MutationAssessor_pred, $MutationAssessor_score,  $MutationAssessor_rank, $CADD_raw, $CADD_rank, $CADD_phred) = (undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($FATHMM_pred, $FATHMM_score, $FATHMM_rank, $MetaSVM_pred, $MetaSVM_score, $MetaSVM_rank, $MetaLR_pred, $MetaLR_score, $MetaLR_rank, $VEST3_score, $VEST3_rank) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
	my ($M_CAP_score, $M_CAP_rank, $M_CAP_pred, $Interpro, $GTEx_V6_gene, $GTEx_V6_tissue, $EIGEN_coding_or_noncoding, $EIGEN_raw, $EIGEN_PC_raw, $GenoCanyon_score, $GenoCanyon_rank) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($GTX_CHIPSEQ, $GTX_CPG, $GTX_DISEASE, $GTX_DNASE, $GTX_DRUG, $GTX_GWAS, $GTX_HGMD, $GTX_HGMD_DISGENE, $GTX_HGMDIMPUTED) = (undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($GTX_MICROSAT, $GTX_MIRNA, $GTX_OMIM, $GTX_ORPHA, $GTX_PATH, $GTX_PGMD, $GTX_PTMS, $GTX_TSS, $GTX_VISTA, $GTX_TRANSFAC_TFBS) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($ENC_HMM, $ENC_TFBS, $SIMPLE_REPEAT, $SEGMENT_DUPS, $REPEAT_MASKER, $GERP_score, $GERP_GT2_WG, $GERP_ELEM_WG) = (undef, undef, undef, undef, undef, undef, undef, undef);
        my ($phyloP100way_vertebrate_score, $phyloP100way_vertebrate_rank, $phyloP20way_mammalian_score, $phyloP20way_mammalian_rank, $phastCons100way_vertebrate_score, $phastCons100way_vertebrate_rank, $phastCons20way_mammalian_score, $phastCons20way_mammalian_rank, $SiPhy_29way_logOdds_score, $SiPhy_29way_logOdds_rank, $WELLDERLY_AF, $WELLDERLY_AN, $REVEL, $INTERVAR, $EBI3222PAK) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($PROVEAN_score, $PROVEAN_pred, $PROVEAN_rank, $DANN_score, $DANN_rank, $FATHMMMKL_coding_score, $FATHMMMKL_coding_rank, $FATHMMMKL_coding_pred, $INTEGRATED_fCons_score, $INTEGRATED_fCons_rank, $INTEGRATED_conf_value, $GERP_RS_dbnsfp, $GERP_rank_dbnsfp) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef); 
        my ($RadialSVM_pred, $LR_pred, $RadialSVM_score, $LR_score, $phyloP46way_placental, $phyloP100way_vertebrate) = (undef, undef, undef, undef, undef, undef);
        my ($RVIS_raw, $RVIS_percent, $ExAC_RVIS_percent, $ExAC_LoF_FDR, $MPC_EXAC ) = (undef, undef, undef, undef, undef);
        my ($UORF_TYPE, $UORF_GENECLASS, $UORF_DISTANCE, $UORF_EFFECT, $UORF_STRENGTH) = (undef, undef, undef, undef, undef);
        my $UV_SIGNATURE = undef;


        foreach my $variantKey (@{$rH_vcfIndex->{$rA_x->[0] . ":" . $rA_x->[1] . "_" . $rA_x->[3] . "/" . $rA_x->[4]}}) {
            $VC                      .= (defined $rH_data->{$variantKey}->{VC})                        ? $rH_data->{$variantKey}->{VC} . "," : "NIL,";
            $VFT                     .= (defined $rH_data->{$variantKey}->{VFT})                       ? join("|", @{$rH_data->{$variantKey}->{VFT}}) . "," : "NIL,";
            $VT                      .= (defined $rH_data->{$variantKey}->{VT})                        ? $rH_data->{$variantKey}->{VT} . "," : "NIL,";
            $DA                      .= (defined $rH_data->{$variantKey}->{DA})                        ? join("|", @{$rH_data->{$variantKey}->{DA}}) . "," : "NIL,";
            $CLINVAR                 .= (defined $rH_data->{$variantKey}->{CLINVAR})                   ? $rH_data->{$variantKey}->{CLINVAR} . "," : "NIL,";
            $COSMIC                  .= (defined $rH_data->{$variantKey}->{COSMIC})                    ? $rH_data->{$variantKey}->{COSMIC} . "," : "NIL,";
            $ICGC21_ID               .= (defined $rH_data->{$variantKey}->{ICGC21_ID})                 ? $rH_data->{$variantKey}->{ICGC21_ID} . "," : "NIL,";
            $ICGC21_Occurrence       .= (defined $rH_data->{$variantKey}->{ICGC21_Occurrence})        ? $rH_data->{$variantKey}->{ICGC21_Occurrence} . "," : "NIL,";
            $NCI60                   .= (defined $rH_data->{$variantKey}->{NCI60})                     ? $rH_data->{$variantKey}->{NCI60} . "," : "NIL,";
            $GTX_GWAS                .= (defined $rH_data->{$variantKey}->{GTX_GWAS})                  ? $rH_data->{$variantKey}->{GTX_GWAS} . "," : "NIL,";
            $GTX_HGMD                .= (defined $rH_data->{$variantKey}->{GTX_HGMD})                  ? $rH_data->{$variantKey}->{GTX_HGMD} . "," : "NIL,";
            $SIFT_pred               .= (defined $rH_data->{$variantKey}->{SIFT_pred})                 ? $rH_data->{$variantKey}->{SIFT_pred} . "," : "NIL,";
            $Polyphen2_HDIV_pred     .= (defined $rH_data->{$variantKey}->{Polyphen2_HDIV_pred})       ? $rH_data->{$variantKey}->{Polyphen2_HDIV_pred} . "," : "NIL,";
            $Polyphen2_HVAR_pred     .= (defined $rH_data->{$variantKey}->{Polyphen2_HVAR_pred})       ? $rH_data->{$variantKey}->{Polyphen2_HVAR_pred} . "," : "NIL,";
            $LRT_pred                .= (defined $rH_data->{$variantKey}->{LRT_pred})                  ? $rH_data->{$variantKey}->{LRT_pred} . "," : "NIL,";
            $MutationTaster_pred     .= (defined $rH_data->{$variantKey}->{MutationTaster_pred})       ? $rH_data->{$variantKey}->{MutationTaster_pred} . "," : "NIL,";
            $MutationAssessor_pred   .= (defined $rH_data->{$variantKey}->{MutationAssessor_pred})     ? $rH_data->{$variantKey}->{MutationAssessor_pred} . "," : "NIL,";
            $FATHMM_pred             .= (defined $rH_data->{$variantKey}->{FATHMM_pred})               ? $rH_data->{$variantKey}->{FATHMM_pred} . "," : "NIL,";
            $RadialSVM_pred          .= (defined $rH_data->{$variantKey}->{RadialSVM_pred})            ? $rH_data->{$variantKey}->{RadialSVM_pred} . "," : "NIL,";
            $LR_pred                 .= (defined $rH_data->{$variantKey}->{LR_pred})                   ? $rH_data->{$variantKey}->{LR_pred} . "," : "NIL,";
            $PROVEAN_pred            .= (defined $rH_data->{$variantKey}->{PROVEAN_pred})              ? $rH_data->{$variantKey}->{PROVEAN_pred} . "," : "NIL,";
            $FATHMMMKL_coding_pred   .= (defined $rH_data->{$variantKey}->{FATHMMMKL_coding_pred})     ? $rH_data->{$variantKey}->{FATHMMMKL_coding_pred} . "," : "NIL,";
            $MetaSVM_pred            .= (defined $rH_data->{$variantKey}->{MetaSVM_pred})              ? $rH_data->{$variantKey}->{MetaSVM_pred} . "," : "NIL,";
            $MetaLR_pred             .= (defined $rH_data->{$variantKey}->{MetaLR_pred})               ? $rH_data->{$variantKey}->{MetaLR_pred} . "," : "NIL,";
            $M_CAP_pred              .= (defined $rH_data->{$variantKey}->{M_CAP_pred})                ? $rH_data->{$variantKey}->{M_CAP_pred} . "," : "NIL,";
            $EIGEN_coding_or_noncoding .= (defined $rH_data->{$variantKey}->{EIGEN_coding_or_noncoding})                ? $rH_data->{$variantKey}->{EIGEN_coding_or_noncoding} . "," : "NIL,";
            $Interpro                .= (defined $rH_data->{$variantKey}->{Interpro})                  ? $rH_data->{$variantKey}->{Interpro} . "," : "NIL,";
            $GTEx_V6_gene            .= (defined $rH_data->{$variantKey}->{GTEx_V6_gene})               ? $rH_data->{$variantKey}->{GTEx_V6_gene} . "," : "NIL,";
            $GTEx_V6_tissue          .= (defined $rH_data->{$variantKey}->{GTEx_V6_tissue})            ? $rH_data->{$variantKey}->{GTEx_V6_tissue} . "," : "NIL,";
            $INTERVAR                .= (defined $rH_data->{$variantKey}->{INTERVAR})                  ? $rH_data->{$variantKey}->{INTERVAR} . "," : "NIL,";
            $oneKG                   .= (defined $rH_data->{$variantKey}->{KG})                        ? $rH_data->{$variantKey}->{KG} . "," : "-1,";
            $AFR                     .= (defined $rH_data->{$variantKey}->{KG_AFR})                    ? $rH_data->{$variantKey}->{KG_AFR} . "," : "-1,";
            $AMR                     .= (defined $rH_data->{$variantKey}->{KG_AMR})                    ? $rH_data->{$variantKey}->{KG_AMR} . "," : "-1,";
            $EUR                     .= (defined $rH_data->{$variantKey}->{KG_EUR})                    ? $rH_data->{$variantKey}->{KG_EUR} . "," : "-1,";
            $EAS                     .= (defined $rH_data->{$variantKey}->{KG_EAS})                    ? $rH_data->{$variantKey}->{KG_EAS} . "," : "-1,";
            $SAS                     .= (defined $rH_data->{$variantKey}->{KG_SAS})                    ? $rH_data->{$variantKey}->{KG_SAS} . "," : "-1,";
            $ExAC_Freq               .= (defined $rH_data->{$variantKey}->{ExAC_Freq})                 ? $rH_data->{$variantKey}->{ExAC_Freq} . "," : "-1,";
            $ExAC_AFR                .= (defined $rH_data->{$variantKey}->{ExAC_AFR})                  ? $rH_data->{$variantKey}->{ExAC_AFR} . "," : "-1,";
            $ExAC_AMR                .= (defined $rH_data->{$variantKey}->{ExAC_AMR})                  ? $rH_data->{$variantKey}->{ExAC_AMR} . "," : "-1,";
            $ExAC_EAS                .= (defined $rH_data->{$variantKey}->{ExAC_EAS})                  ? $rH_data->{$variantKey}->{ExAC_EAS} . "," : "-1,";
            $ExAC_FIN                .= (defined $rH_data->{$variantKey}->{ExAC_FIN})                  ? $rH_data->{$variantKey}->{ExAC_FIN} . "," : "-1,";
            $ExAC_NFE                .= (defined $rH_data->{$variantKey}->{ExAC_NFE})                  ? $rH_data->{$variantKey}->{ExAC_NFE} . "," : "-1,";
            $ExAC_OTH                .= (defined $rH_data->{$variantKey}->{ExAC_OTH})                  ? $rH_data->{$variantKey}->{ExAC_OTH} . "," : "-1,";
            $ExAC_SAS                .= (defined $rH_data->{$variantKey}->{ExAC_SAS})                  ? $rH_data->{$variantKey}->{ExAC_SAS} . "," : "-1,";
            $GME_AF                  .= (defined $rH_data->{$variantKey}->{GME_AF})                    ? $rH_data->{$variantKey}->{GME_AF} . "," : "-1,";
            $GME_NWA                 .= (defined $rH_data->{$variantKey}->{GME_NWA})                   ? $rH_data->{$variantKey}->{GME_NWA} . "," : "-1,";
            $GME_NEA                 .= (defined $rH_data->{$variantKey}->{GME_NEA})                   ? $rH_data->{$variantKey}->{GME_NEA} . "," : "-1,"; 
            $GME_AP                  .= (defined $rH_data->{$variantKey}->{GME_AP})                    ? $rH_data->{$variantKey}->{GME_AP} . "," : "-1,";
            $GME_Israel              .= (defined $rH_data->{$variantKey}->{GME_Israel})                ? $rH_data->{$variantKey}->{GME_Israel} . "," : "-1,"; 
            $GME_SD                  .= (defined $rH_data->{$variantKey}->{GME_SD})                    ? $rH_data->{$variantKey}->{GME_SD} . "," : "-1,"; 
            $GME_TP                  .= (defined $rH_data->{$variantKey}->{GME_TP})                    ? $rH_data->{$variantKey}->{GME_TP} . "," : "-1,";
            $GME_CA                  .= (defined $rH_data->{$variantKey}->{GME_CA})                    ? $rH_data->{$variantKey}->{GME_CA} . "," : "-1,";
            $gnomAD_exome_ALL        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_ALL})          ? $rH_data->{$variantKey}->{gnomAD_exome_ALL} . "," : "-1,";
            $gnomAD_exome_AFR        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_AFR})          ? $rH_data->{$variantKey}->{gnomAD_exome_AFR} . "," : "-1,";
            $gnomAD_exome_AMR        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_AMR})          ? $rH_data->{$variantKey}->{gnomAD_exome_AMR} . "," : "-1,";
            $gnomAD_exome_ASJ        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_ASJ})          ? $rH_data->{$variantKey}->{gnomAD_exome_ASJ} . "," : "-1,";
            $gnomAD_exome_EAS        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_EAS})          ? $rH_data->{$variantKey}->{gnomAD_exome_EAS} . "," : "-1,";
            $gnomAD_exome_FIN        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_FIN})          ? $rH_data->{$variantKey}->{gnomAD_exome_FIN} . "," : "-1,";
            $gnomAD_exome_NFE        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_NFE})          ? $rH_data->{$variantKey}->{gnomAD_exome_NFE} . "," : "-1,";
            $gnomAD_exome_OTH        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_OTH})          ? $rH_data->{$variantKey}->{gnomAD_exome_OTH} . "," : "-1,";
            $gnomAD_exome_SAS        .= (defined $rH_data->{$variantKey}->{gnomAD_exome_SAS})          ? $rH_data->{$variantKey}->{gnomAD_exome_SAS} . "," : "-1,";
            $gnomAD_genome_ALL       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_ALL})         ? $rH_data->{$variantKey}->{gnomAD_genome_ALL} . "," : "-1,";
            $gnomAD_genome_AFR       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_AFR})         ? $rH_data->{$variantKey}->{gnomAD_genome_AFR} . "," : "-1,";
            $gnomAD_genome_AMR       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_AMR})         ? $rH_data->{$variantKey}->{gnomAD_genome_AMR} . "," : "-1,";
            $gnomAD_genome_ASJ       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_ASJ})         ? $rH_data->{$variantKey}->{gnomAD_genome_ASJ} . "," : "-1,";
            $gnomAD_genome_EAS       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_EAS})         ? $rH_data->{$variantKey}->{gnomAD_genome_EAS} . "," : "-1,";
            $gnomAD_genome_FIN       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_FIN})         ? $rH_data->{$variantKey}->{gnomAD_genome_FIN} . "," : "-1,";
            $gnomAD_genome_NFE       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_NFE})         ? $rH_data->{$variantKey}->{gnomAD_genome_NFE} . "," : "-1,";
            $gnomAD_genome_OTH       .= (defined $rH_data->{$variantKey}->{gnomAD_genome_OTH})         ? $rH_data->{$variantKey}->{gnomAD_genome_OTH} . "," : "-1,";
            $GoNL                    .= (defined $rH_data->{$variantKey}->{GoNL})                      ? $rH_data->{$variantKey}->{GoNL} . "," : "-1,";
            $WELLDERLY_AF            .= (defined $rH_data->{$variantKey}->{WELLDERLY_AF})              ? $rH_data->{$variantKey}->{WELLDERLY_AF} . "," : "-1,";
            $WELLDERLY_AN            .= (defined $rH_data->{$variantKey}->{WELLDERLY_AN})              ? $rH_data->{$variantKey}->{WELLDERLY_AN} . "," : "-1,";
            $CG46                    .= (defined $rH_data->{$variantKey}->{CG46})                      ? $rH_data->{$variantKey}->{CG46} . "," : "-1,";
            $ESP6500                 .= (defined $rH_data->{$variantKey}->{ESP6500})                   ? $rH_data->{$variantKey}->{ESP6500} . "," : "-1,";
            $ESP6500AA               .= (defined $rH_data->{$variantKey}->{ESP6500AA})                 ? $rH_data->{$variantKey}->{ESP6500AA} . "," : "-1,";
            $ESP6500EA               .= (defined $rH_data->{$variantKey}->{ESP6500EA})                 ? $rH_data->{$variantKey}->{ESP6500EA} . "," : "-1,";
            $SIFT_score              .= (defined $rH_data->{$variantKey}->{SIFT_score})                ? $rH_data->{$variantKey}->{SIFT_score} . "," : "-1,";
            $SIFT_rank               .= (defined $rH_data->{$variantKey}->{SIFT_rank})                 ? $rH_data->{$variantKey}->{SIFT_rank} . "," : "-1,";
            $Polyphen2_HDIV_score    .= (defined $rH_data->{$variantKey}->{Polyphen2_HDIV_score})      ? $rH_data->{$variantKey}->{Polyphen2_HDIV_score} . "," : "-1,";
            $Polyphen2_HDIV_rank     .= (defined $rH_data->{$variantKey}->{Polyphen2_HDIV_rank})       ? $rH_data->{$variantKey}->{Polyphen2_HDIV_rank} . "," : "-1,";
            $Polyphen2_HVAR_score    .= (defined $rH_data->{$variantKey}->{Polyphen2_HVAR_score})      ? $rH_data->{$variantKey}->{Polyphen2_HVAR_score} . "," : "-1,";
            $Polyphen2_HVAR_rank     .= (defined $rH_data->{$variantKey}->{Polyphen2_HVAR_rank})       ? $rH_data->{$variantKey}->{Polyphen2_HVAR_rank} . "," : "-1,";
            $LRT_score               .= (defined $rH_data->{$variantKey}->{LRT_score})                 ? $rH_data->{$variantKey}->{LRT_score} . "," : "-1,";
            $LRT_rank               .= (defined $rH_data->{$variantKey}->{LRT_rank})                 ? $rH_data->{$variantKey}->{LRT_rank} . "," : "-1,";
            $MutationTaster_score    .= (defined $rH_data->{$variantKey}->{MutationTaster_score})      ? $rH_data->{$variantKey}->{MutationTaster_score} . "," : "-1,";
            $MutationTaster_rank    .= (defined $rH_data->{$variantKey}->{MutationTaster_rank})      ? $rH_data->{$variantKey}->{MutationTaster_rank} . "," : "-1,";
            $MutationAssessor_score  .= (defined $rH_data->{$variantKey}->{MutationAssessor_score})    ? $rH_data->{$variantKey}->{MutationAssessor_score} . "," : "-1,";
            $MutationAssessor_rank  .= (defined $rH_data->{$variantKey}->{MutationAssessor_rank})    ? $rH_data->{$variantKey}->{MutationAssessor_rank} . "," : "-1,";
            $FATHMM_score            .= (defined $rH_data->{$variantKey}->{FATHMM_score})              ? $rH_data->{$variantKey}->{FATHMM_score} . "," : "-1,";
            $FATHMM_rank            .= (defined $rH_data->{$variantKey}->{FATHMM_rank})              ? $rH_data->{$variantKey}->{FATHMM_rank} . "," : "-1,";
            $PROVEAN_score           .= (defined $rH_data->{$variantKey}->{PROVEAN_score})             ? $rH_data->{$variantKey}->{PROVEAN_score} . "," : "-1,";
            $PROVEAN_rank           .= (defined $rH_data->{$variantKey}->{PROVEAN_rank})             ? $rH_data->{$variantKey}->{PROVEAN_rank} . "," : "-1,";
            $VEST3_score             .= (defined $rH_data->{$variantKey}->{VEST3_score})               ? $rH_data->{$variantKey}->{VEST3_score} . "," : "-1,";
            $VEST3_rank             .= (defined $rH_data->{$variantKey}->{VEST3_rank})               ? $rH_data->{$variantKey}->{VEST3_rank} . "," : "-1,";
            $MetaSVM_score           .= (defined $rH_data->{$variantKey}->{MetaSVM_score})             ? $rH_data->{$variantKey}->{MetaSVM_score} . "," : "-1,";
            $MetaSVM_rank           .= (defined $rH_data->{$variantKey}->{MetaSVM_rank})             ? $rH_data->{$variantKey}->{MetaSVM_rank} . "," : "-1,";
            $MetaLR_score            .= (defined $rH_data->{$variantKey}->{MetaLR_score})              ? $rH_data->{$variantKey}->{MetaLR_score} . "," : "-1,";
            $MetaLR_rank            .= (defined $rH_data->{$variantKey}->{MetaLR_rank})              ? $rH_data->{$variantKey}->{MetaLR_rank} . "," : "-1,";
            $M_CAP_score            .= (defined $rH_data->{$variantKey}->{M_CAP_score})              ? $rH_data->{$variantKey}->{M_CAP_score} . "," : "-1,";
            $M_CAP_rank            .= (defined $rH_data->{$variantKey}->{M_CAP_rank})              ? $rH_data->{$variantKey}->{M_CAP_rank} . "," : "-1,";
            $CADD_raw                .= (defined $rH_data->{$variantKey}->{CADD_raw})                  ? $rH_data->{$variantKey}->{CADD_raw} . "," : "-1,";
            $CADD_phred              .= (defined $rH_data->{$variantKey}->{CADD_phred})                ? $rH_data->{$variantKey}->{CADD_phred} . "," : "-1,";
            $CADD_rank              .= (defined $rH_data->{$variantKey}->{CADD_rank})                ? $rH_data->{$variantKey}->{CADD_rank} . "," : "-1,";
            $DANN_score            .= (defined $rH_data->{$variantKey}->{DANN_score})              ? $rH_data->{$variantKey}->{DANN_score} . "," : "-1,";
            $DANN_rank            .= (defined $rH_data->{$variantKey}->{DANN_rank})              ? $rH_data->{$variantKey}->{DANN_rank} . "," : "-1,";
            $FATHMMMKL_coding_score  .= (defined $rH_data->{$variantKey}->{FATHMMMKL_coding_score})    ? $rH_data->{$variantKey}->{FATHMMMKL_coding_score} . "," : "-1,";
            $FATHMMMKL_coding_rank  .= (defined $rH_data->{$variantKey}->{FATHMMMKL_coding_rank})    ? $rH_data->{$variantKey}->{FATHMMMKL_coding_rank} . "," : "-1,";
            $EIGEN_raw                .= (defined $rH_data->{$variantKey}->{EIGEN_raw})                  ? $rH_data->{$variantKey}->{EIGEN_raw} . "," : "-1,";
            $EIGEN_PC_raw                .= (defined $rH_data->{$variantKey}->{EIGEN_PC_raw})                  ? $rH_data->{$variantKey}->{EIGEN_PC_raw} . "," : "-1,";
            $GenoCanyon_score            .= (defined $rH_data->{$variantKey}->{GenoCanyon_score})              ? $rH_data->{$variantKey}->{GenoCanyon_score} . "," : "-1,";
            $GenoCanyon_rank            .= (defined $rH_data->{$variantKey}->{GenoCanyon_rank})              ? $rH_data->{$variantKey}->{GenoCanyon_rank} . "," : "-1,";
            $RadialSVM_score         .= (defined $rH_data->{$variantKey}->{RadialSVM_score})           ? $rH_data->{$variantKey}->{RadialSVM_score} . "," : "-1,";
            $LR_score                .= (defined $rH_data->{$variantKey}->{LR_score})                  ? $rH_data->{$variantKey}->{LR_score} . "," : "-1,";
            $INTEGRATED_fCons_score  .= (defined $rH_data->{$variantKey}->{INTEGRATED_fCons_score})    ? $rH_data->{$variantKey}->{INTEGRATED_fCons_score} . "," : "-1,";
            $INTEGRATED_fCons_rank  .= (defined $rH_data->{$variantKey}->{INTEGRATED_fCons_rank})    ? $rH_data->{$variantKey}->{INTEGRATED_fCons_rank} . "," : "-1,";
            $INTEGRATED_conf_value   .= (defined $rH_data->{$variantKey}->{INTEGRATED_conf_value})     ? $rH_data->{$variantKey}->{INTEGRATED_conf_value} . "," : "-1,";
            $GERP_RS_dbnsfp                 .= (defined $rH_data->{$variantKey}->{GERP_RS_dbnsfp})               ? $rH_data->{$variantKey}->{GERP_RS_dbnsfp} . "," : "-1,";
            $GERP_rank_dbnsfp                 .= (defined $rH_data->{$variantKey}->{GERP_rank_dbnsfp})               ? $rH_data->{$variantKey}->{GERP_rank_dbnsfp} . "," : "-1,";
            $GERP_GT2_WG                  .= (defined $rH_data->{$variantKey}->{GERP_GT2_WG})               ? $rH_data->{$variantKey}->{GERP_GT2_WG} . "," : "-1,";
            $phyloP100way_vertebrate_score   .= (defined $rH_data->{$variantKey}->{phyloP100way_vertebrate_score})     ? $rH_data->{$variantKey}->{phyloP100way_vertebrate_score} . "," : "-1,";
            $phyloP100way_vertebrate_rank   .= (defined $rH_data->{$variantKey}->{phyloP100way_vertebrate_rank})     ? $rH_data->{$variantKey}->{phyloP100way_vertebrate_rank} . "," : "-1,";
            $phyloP20way_mammalian_score   .= (defined $rH_data->{$variantKey}->{phyloP20way_mammalian_score})     ? $rH_data->{$variantKey}->{phyloP20way_mammalian_score} . "," : "-1,";
            $phyloP20way_mammalian_rank   .= (defined $rH_data->{$variantKey}->{phyloP20way_mammalian_rank})     ? $rH_data->{$variantKey}->{phyloP20way_mammalian_rank} . "," : "-1,";
            $phastCons100way_vertebrate_score .= (defined $rH_data->{$variantKey}->{phastCons100way_vertebrate_score}) ? $rH_data->{$variantKey}->{phastCons100way_vertebrate_score} . "," : "-1,";
            $phastCons100way_vertebrate_rank .= (defined $rH_data->{$variantKey}->{phastCons100way_vertebrate_rank}) ? $rH_data->{$variantKey}->{phastCons100way_vertebrate_rank} . "," : "-1,";
            $phastCons20way_mammalian_score .= (defined $rH_data->{$variantKey}->{phastCons20way_mammalian_score}) ? $rH_data->{$variantKey}->{phastCons20way_mammalian_score} . "," : "-1,";
            $phastCons20way_mammalian_rank .= (defined $rH_data->{$variantKey}->{phastCons20way_mammalian_rank}) ? $rH_data->{$variantKey}->{phastCons20way_mammalian_rank} . "," : "-1,";
            $SiPhy_29way_logOdds_score    .= (defined $rH_data->{$variantKey}->{SiPhy_29way_logOdds_score})       ? $rH_data->{$variantKey}->{SiPhy_29way_logOdds_score} . "," : "-1,";
            $SiPhy_29way_logOdds_rank    .= (defined $rH_data->{$variantKey}->{SiPhy_29way_logOdds_rank})       ? $rH_data->{$variantKey}->{SiPhy_29way_logOdds_rank} . "," : "-1,";
            $REVEL                   .= (defined $rH_data->{$variantKey}->{REVEL})                     ? $rH_data->{$variantKey}->{REVEL}. "," : "-1,";
	    $WGSCADD                 .= (defined $rH_data->{$variantKey}->{WGSCADD})                   ? $rH_data->{$variantKey}->{WGSCADD} . "," : "-1,";
	    $WGSCADD_Phred           .= (defined $rH_data->{$variantKey}->{WGSCADD_Phred})             ? $rH_data->{$variantKey}->{WGSCADD_Phred} . "," : "-1,";
            $WGSDANN                 .= (defined $rH_data->{$variantKey}->{WGSDANN})                   ? $rH_data->{$variantKey}->{WGSDANN} . "," : "-1,";
            $dbscSNV_ADA_SCORE       .= (defined $rH_data->{$variantKey}->{dbscSNV_ADA_SCORE})         ? $rH_data->{$variantKey}->{dbscSNV_ADA_SCORE} . "," : "-1,";
            $dbscSNV_RF_SCORE        .= (defined $rH_data->{$variantKey}->{dbscSNV_RF_SCORE})          ? $rH_data->{$variantKey}->{dbscSNV_RF_SCORE} . "," : "-1,";
            $HRCR1_AF                .= (defined $rH_data->{$variantKey}->{HRCR1_AF})                  ? $rH_data->{$variantKey}->{HRCR1_AF} . "," : "-1,";
            $HRCR1_AC                .= (defined $rH_data->{$variantKey}->{HRCR1_AC})                  ? $rH_data->{$variantKey}->{HRCR1_AC} . "," : "-1,";
            $HRCR1_AN                .= (defined $rH_data->{$variantKey}->{HRCR1_AN})                  ? $rH_data->{$variantKey}->{HRCR1_AN} . "," : "-1,";
            $HRCR1_nonKG_AF          .= (defined $rH_data->{$variantKey}->{HRCR1_nonKG_AF})            ? $rH_data->{$variantKey}->{HRCR1_nonKG_AF} . "," : "-1,";
            $HRCR1_nonKG_AC          .= (defined $rH_data->{$variantKey}->{HRCR1_nonKG_AC})            ? $rH_data->{$variantKey}->{HRCR1_nonKG_AC} . "," : "-1,";
            $HRCR1_nonKG_AN          .= (defined $rH_data->{$variantKey}->{HRCR1_nonKG_AN})            ? $rH_data->{$variantKey}->{HRCR1_nonKG_AN} . "," : "-1,";
            $WGSFATHMM_coding        .= (defined $rH_data->{$variantKey}->{WGSFATHMM_coding})          ? $rH_data->{$variantKey}->{WGSFATHMM_coding} . "," : "-1,";
            $WGSFATHMM_noncoding     .= (defined $rH_data->{$variantKey}->{WGSFATHMM_noncoding})       ? $rH_data->{$variantKey}->{WGSFATHMM_noncoding} . "," : "-1,";
            $GWAVA_region_score      .= (defined $rH_data->{$variantKey}->{GWAVA_region_score})        ? $rH_data->{$variantKey}->{GWAVA_region_score} . "," : "-1,";
            $GWAVA_tss_score         .= (defined $rH_data->{$variantKey}->{GWAVA_tss_score})           ? $rH_data->{$variantKey}->{GWAVA_tss_score} . "," : "-1,";
	    $ExAC_ALL_nonpsych       .= (defined $rH_data->{$variantKey}->{ExAC_ALL_nonpsych})         ? $rH_data->{$variantKey}->{ExAC_ALL_nonpsych} . "," : "-1,";
	    $ExAC_ALL_nontcga        .= (defined $rH_data->{$variantKey}->{ExAC_ALL_nontcga})          ? $rH_data->{$variantKey}->{ExAC_ALL_nontcga} . "," : "-1,";
            $KAVIAR_AF               .= (defined $rH_data->{$variantKey}->{KAVIAR_AF})                 ? $rH_data->{$variantKey}->{KAVIAR_AF} . "," : "-1,";
            $KAVIAR_AC               .= (defined $rH_data->{$variantKey}->{KAVIAR_AC})                 ? $rH_data->{$variantKey}->{KAVIAR_AC} . "," : "-1,";
            $KAVIAR_AN               .= (defined $rH_data->{$variantKey}->{KAVIAR_AN})                 ? $rH_data->{$variantKey}->{KAVIAR_AN} . "," : "-1,";
            $WGSEIGEN	             .= (defined $rH_data->{$variantKey}->{WGSEIGEN}) 		       ? $rH_data->{$variantKey}->{WGSEIGEN} . "," : "-1,";
            $EBI3222PAK              .= (defined $rH_data->{$variantKey}->{EBI3222PAK})                ? $rH_data->{$variantKey}->{EBI3222PAK} . "," : "-1";
            $MPC_EXAC                .= (defined $rH_data->{$variantKey}->{MPC_EXAC})                  ? $rH_data->{$variantKey}->{MPC_EXAC} . "," : "-1";
            $UORF_TYPE               .= (defined $rH_data->{$variantKey}->{uORF_type})                 ? $rH_data->{$variantKey}->{uORF_type} . "," : "-1";
            $UORF_GENECLASS              .= (defined $rH_data->{$variantKey}->{uORF_geneClass})                ? $rH_data->{$variantKey}->{uORF_geneClass} . "," : "-1";
            $UORF_DISTANCE           .= (defined $rH_data->{$variantKey}->{uORF_distance})             ? $rH_data->{$variantKey}->{uORF_distance} . "," : "-1";
            $UORF_EFFECT             .= (defined $rH_data->{$variantKey}->{uORF_effect})               ? $rH_data->{$variantKey}->{uORF_effect} . "," : "-1";
            $UORF_STRENGTH           .= (defined $rH_data->{$variantKey}->{uORF_kozakStrength})             ? $rH_data->{$variantKey}->{uORF_kozakStrength} . "," : "-1";
            $UV_SIGNATURE            .= (defined $rH_data->{$variantKey}->{UV_signature})                 ? $rH_data->{$variantKey}->{UV_signature} . "," : "-1";


            ((defined $rH_data->{$variantKey}->{Gene}) && (!defined($Gene)))                                       and $Gene                    .= join("|", @{$rH_data->{$variantKey}->{Gene}}) . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_132}) && (!defined($dbSNP_132)))                             and $dbSNP_132               .= $rH_data->{$variantKey}->{dbSNP_132} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_138}) && (!defined($dbSNP_138)))                             and $dbSNP_138               .= $rH_data->{$variantKey}->{dbSNP_138} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_138Common}) && (!defined($dbSNP_138Common)))                 and $dbSNP_138Common         .= $rH_data->{$variantKey}->{dbSNP_138Common} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_138CommonAlCount}) && (!defined($dbSNP_138CommonAlCount)))   and $dbSNP_138CommonAlCount  .= $rH_data->{$variantKey}->{dbSNP_138CommonAlCount} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_138CommonAlFreq}) && (!defined($dbSNP_138CommonAlFreq)))     and $dbSNP_138CommonAlFreq   .= $rH_data->{$variantKey}->{dbSNP_138CommonAlFreq} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_138Flagged}) && (!defined($dbSNP_138Flagged)))               and $dbSNP_138Flagged        .= $rH_data->{$variantKey}->{dbSNP_138Flagged} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_138FlaggedAlCount}) && (!defined($dbSNP_138FlaggedAlCount))) and $dbSNP_138FlaggedAlCount .= $rH_data->{$variantKey}->{dbSNP_138FlaggedAlCount} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_138FlaggedAlFreq}) && (!defined($dbSNP_138FlaggedAlFreq)))   and $dbSNP_138FlaggedAlFreq  .= $rH_data->{$variantKey}->{dbSNP_138FlaggedAlFreq} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_142}) && (!defined($dbSNP_142)))                             and $dbSNP_142               .= $rH_data->{$variantKey}->{dbSNP_142} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_142Common}) && (!defined($dbSNP_142Common)))                 and $dbSNP_142Common         .= $rH_data->{$variantKey}->{dbSNP_142Common} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_142CommonAlCount}) && (!defined($dbSNP_142CommonAlCount)))   and $dbSNP_142CommonAlCount  .= $rH_data->{$variantKey}->{dbSNP_142CommonAlCount} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_142CommonAlFreq}) && (!defined($dbSNP_142CommonAlFreq)))     and $dbSNP_142CommonAlFreq   .= $rH_data->{$variantKey}->{dbSNP_142CommonAlFreq} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_142Flagged}) && (!defined($dbSNP_142Flagged)))               and $dbSNP_142Flagged        .= $rH_data->{$variantKey}->{dbSNP_142Flagged} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_142FlaggedAlCount}) && (!defined($dbSNP_142FlaggedAlCount))) and $dbSNP_142FlaggedAlCount .= $rH_data->{$variantKey}->{dbSNP_142FlaggedAlCount} . ",";
            ((defined $rH_data->{$variantKey}->{dbSNP_142FlaggedAlFreq}) && (!defined($dbSNP_142FlaggedAlFreq)))   and $dbSNP_142FlaggedAlFreq  .= $rH_data->{$variantKey}->{dbSNP_142FlaggedAlFreq} . ",";
            ((defined $rH_data->{$variantKey}->{GERP_ELEM_WG}) && (!defined($GERP_ELEM_WG)))                           and $GERP_ELEM_WG                 = $rH_data->{$variantKey}->{GERP_ELEM_WG} . ",";
            ((defined $rH_data->{$variantKey}->{MIRNA}) && (!defined($MIRNA)))                                     and $MIRNA                   .= $rH_data->{$variantKey}->{MIRNA} . ",";
            ((defined $rH_data->{$variantKey}->{MIRNATARGET}) && (!defined($MIRNATARGET)))                         and $MIRNATARGET             .= $rH_data->{$variantKey}->{MIRNATARGET} . ",";
            ((defined $rH_data->{$variantKey}->{RNASEQ_CUFF}) && (!defined($RNASEQ_CUFF)))                         and $RNASEQ_CUFF             .= $rH_data->{$variantKey}->{RNASEQ_CUFF} . ",";
            ((defined $rH_data->{$variantKey}->{RNASEQ_SCRI}) && (!defined($RNASEQ_SCRI)))                         and $RNASEQ_SCRI             .= $rH_data->{$variantKey}->{RNASEQ_SCRI} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_CHIPSEQ}) && (!defined($GTX_CHIPSEQ)))                         and $GTX_CHIPSEQ             .= $rH_data->{$variantKey}->{GTX_CHIPSEQ} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_CPG}) && (!defined($GTX_CPG)))                                 and $GTX_CPG                 .= $rH_data->{$variantKey}->{GTX_CPG} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_DISEASE}) && (!defined($GTX_DISEASE)))                         and $GTX_DISEASE             .= $rH_data->{$variantKey}->{GTX_DISEASE} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_DNASE}) && (!defined($GTX_DNASE)))                             and $GTX_DNASE               .= $rH_data->{$variantKey}->{GTX_DNASE} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_DRUG}) && (!defined($GTX_DRUG)))                               and $GTX_DRUG                .= $rH_data->{$variantKey}->{GTX_DRUG} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_HGMD_DISGENE}) && (!defined($GTX_HGMD_DISGENE)))               and $GTX_HGMD_DISGENE        .= $rH_data->{$variantKey}->{GTX_HGMD_DISGENE} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_HGMDIMPUTED}) && (!defined($GTX_HGMDIMPUTED)))                 and $GTX_HGMDIMPUTED         .= $rH_data->{$variantKey}->{GTX_HGMDIMPUTED} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_MICROSAT}) && (!defined($GTX_MICROSAT)))                       and $GTX_MICROSAT            .= $rH_data->{$variantKey}->{GTX_MICROSAT} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_MIRNA}) && (!defined($GTX_MIRNA)))                             and $GTX_MIRNA               .= $rH_data->{$variantKey}->{GTX_MIRNA} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_OMIM}) && (!defined($GTX_OMIM)))                               and $GTX_OMIM                .= $rH_data->{$variantKey}->{GTX_OMIM} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_ORPHA}) && (!defined($GTX_ORPHA)))                             and $GTX_ORPHA               .= $rH_data->{$variantKey}->{GTX_ORPHA} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_PATH}) && (!defined($GTX_PATH)))                               and $GTX_PATH                .= $rH_data->{$variantKey}->{GTX_PATH} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_PGMD}) && (!defined($GTX_PGMD)))                               and $GTX_PGMD                 .= $rH_data->{$variantKey}->{GTX_PGMD} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_PTMS}) && (!defined($GTX_PTMS)))                               and $GTX_PTMS                .= $rH_data->{$variantKey}->{GTX_PTMS} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_TRANSFAC_TFBS}) && (!defined($GTX_TRANSFAC_TFBS)))             and $GTX_TRANSFAC_TFBS       .= $rH_data->{$variantKey}->{GTX_TRANSFAC_TFBS} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_TSS}) && (!defined($GTX_TSS)))                                 and $GTX_TSS                 .= $rH_data->{$variantKey}->{GTX_TSS} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_VISTA}) && (!defined($GTX_VISTA)))                             and $GTX_VISTA               .= $rH_data->{$variantKey}->{GTX_VISTA} . ",";
            ((defined $rH_data->{$variantKey}->{ENC_HMM}) && (!defined($ENC_HMM)))                                 and $ENC_HMM                 .= $rH_data->{$variantKey}->{ENC_HMM} . ",";
            ((defined $rH_data->{$variantKey}->{ENC_TFBS}) && (!defined($ENC_TFBS)))                               and $ENC_TFBS                .= $rH_data->{$variantKey}->{ENC_TFBS} . ",";
            ((defined $rH_data->{$variantKey}->{SIMPLE_REPEAT}) && (!defined($SIMPLE_REPEAT)))                     and $SIMPLE_REPEAT           .= $rH_data->{$variantKey}->{SIMPLE_REPEAT} . ",";
            ((defined $rH_data->{$variantKey}->{SEGMENT_DUPS}) && (!defined($SEGMENT_DUPS)))                       and $SEGMENT_DUPS            .= $rH_data->{$variantKey}->{SEGMENT_DUPS} . ",";
            ((defined $rH_data->{$variantKey}->{REPEAT_MASKER}) && (!defined($REPEAT_MASKER)))                     and $REPEAT_MASKER           .= $rH_data->{$variantKey}->{REPEAT_MASKER} . ",";
            ((defined $rH_data->{$variantKey}->{RVIS_raw}) && (!defined($RVIS_raw)))                               and $RVIS_raw            .= $rH_data->{$variantKey}->{RVIS_raw} . ",";
            ((defined $rH_data->{$variantKey}->{RVIS_percent}) && (!defined($RVIS_percent)))                       and $RVIS_percent            .= $rH_data->{$variantKey}->{RVIS_percent} . ",";
            ((defined $rH_data->{$variantKey}->{ExAC_RVIS_percent}) && (!defined($ExAC_RVIS_percent)))             and $ExAC_RVIS_percent            .= $rH_data->{$variantKey}->{ExAC_RVIS_percent} . ",";
            ((defined $rH_data->{$variantKey}->{ExAC_LoF_FDR}) && (!defined($ExAC_LoF_FDR)))                       and $ExAC_LoF_FDR            .= $rH_data->{$variantKey}->{ExAC_LoF_FDR} . ",";
            ((defined $rH_data->{$variantKey}->{MPC_EXAC}) && (!defined($MPC_EXAC)))                               and $MPC_EXAC                .= $rH_data->{$variantKey}->{MPC_EXAC} . ",";
            ((defined $rH_data->{$variantKey}->{uORF_type}) && (!defined($UORF_TYPE)))                             and $UORF_TYPE                .= $rH_data->{$variantKey}->{uORF_type} . ",";
            ((defined $rH_data->{$variantKey}->{uORF_geneClass}) && (!defined($UORF_GENECLASS)))                           and $UORF_GENECLASS                .= $rH_data->{$variantKey}->{uORF_geneClass} . ",";
            ((defined $rH_data->{$variantKey}->{uORF_distance}) && (!defined($UORF_DISTANCE)))                     and $UORF_DISTANCE                .= $rH_data->{$variantKey}->{uORF_distance} . ",";
            ((defined $rH_data->{$variantKey}->{uORF_effect}) && (!defined($UORF_EFFECT)))                         and $UORF_EFFECT                .= $rH_data->{$variantKey}->{uORF_effect} . ",";
            ((defined $rH_data->{$variantKey}->{uORF_kozakStrength}) && (!defined($UORF_STRENGTH)))                     and $UORF_STRENGTH                .= $rH_data->{$variantKey}->{uORF_kozakStrength} . ",";
            ((defined $rH_data->{$variantKey}->{UV_signature}) && (!defined($UV_SIGNATURE)))                             and $UV_SIGNATURE  .= $rH_data->{$variantKey}->{UV_signature} . ",";

        }

        # Remove all trailing comma and spaces from the annotation strings
        foreach my $annot ( $VC, $VFT, $VT, $Gene, $DA, $dbSNP_132, $dbSNP_138, $dbSNP_138Common, $dbSNP_138CommonAlCount, $dbSNP_138CommonAlFreq, $dbSNP_138Flagged,
                            $dbSNP_138FlaggedAlCount, $dbSNP_138FlaggedAlFreq,$dbSNP_142, $dbSNP_142Common, $dbSNP_142CommonAlCount, $dbSNP_142CommonAlFreq, $dbSNP_142Flagged,
                            $dbSNP_142FlaggedAlCount, $dbSNP_142FlaggedAlFreq,$oneKG, $CG46, $SIFT_score, $SIFT_rank, $SIFT_pred, $Polyphen2_HDIV_score, $Polyphen2_HDIV_rank, $Polyphen2_HDIV_pred,
                            $Polyphen2_HVAR_score, $Polyphen2_HVAR_rank, $Polyphen2_HVAR_pred, $LRT_score, $LRT_rank, $LRT_pred, $MutationTaster_score, $MutationTaster_rank, $MutationTaster_pred, 
                            $MutationAssessor_score, $MutationAssessor_rank, $M_CAP_score, $M_CAP_rank, $M_CAP_pred,
                            $MutationAssessor_pred, $FATHMM_score, $FATHMM_rank, $FATHMM_pred, $RadialSVM_score, $RadialSVM_pred, $LR_score, $LR_pred, $VEST3_score, $VEST3_rank, 
                            $CADD_raw, $CADD_rank, $CADD_phred, $DANN_score, $DANN_rank, $EIGEN_coding_or_noncoding, $EIGEN_raw, $EIGEN_PC_raw, $GenoCanyon_score, $GenoCanyon_rank,
                            $PROVEAN_score, $PROVEAN_rank, $PROVEAN_pred, $FATHMMMKL_coding_score, $FATHMMMKL_coding_rank, $FATHMMMKL_coding_pred, 
                            $INTEGRATED_fCons_score, $INTEGRATED_fCons_rank, $INTEGRATED_conf_value, $INTERVAR, $EBI3222PAK,
                            $MetaSVM_score, $MetaSVM_rank, $MetaSVM_pred, $MetaLR_score, $MetaLR_pred, $MetaLR_rank, $GERP_RS_dbnsfp, $GERP_rank_dbnsfp, 
                            $GERP_score, $GERP_GT2_WG, $GERP_ELEM_WG, $REVEL, $Interpro, $GTEx_V6_gene, $GTEx_V6_tissue,
                            $phyloP100way_vertebrate_score, $phyloP100way_vertebrate_rank, $phyloP20way_mammalian_score, $phyloP20way_mammalian_rank, 
                            $phastCons100way_vertebrate_score, $phastCons100way_vertebrate_rank, $phastCons20way_mammalian_score, $phastCons20way_mammalian_rank,
                            $SiPhy_29way_logOdds_score, $SiPhy_29way_logOdds_rank, $ESP6500, $ESP6500AA, $ESP6500EA,
                            $gnomAD_exome_ALL, $gnomAD_exome_AFR, $gnomAD_exome_AMR, $gnomAD_exome_ASJ, $gnomAD_exome_EAS, $gnomAD_exome_FIN, $gnomAD_exome_NFE, $gnomAD_exome_OTH, $gnomAD_exome_SAS,
                            $gnomAD_genome_ALL, $gnomAD_genome_AFR, $gnomAD_genome_AMR, $gnomAD_genome_ASJ, $gnomAD_genome_EAS, $gnomAD_genome_FIN, $gnomAD_genome_NFE, $gnomAD_genome_OTH, 
                            $MIRNA, $MIRNATARGET, $AFR, $AMR, $EUR, $EAS, $SAS, $GoNL, $RNASEQ_CUFF, $RNASEQ_SCRI, $GTX_CHIPSEQ, $CLINVAR, $COSMIC, $ICGC21_ID, $ICGC21_Occurrence, $GTX_CPG,
                            $GTX_TSS, $GTX_VISTA, $GTX_DISEASE, $GTX_DNASE, $GTX_DRUG, $GTX_GWAS, $GTX_HGMD, $GTX_HGMD_DISGENE, $GTX_HGMDIMPUTED, $GTX_MICROSAT, $GTX_MIRNA, 
                            $GTX_OMIM, $GTX_ORPHA, $GTX_PATH, $GTX_PGMD, $GTX_PTMS,
                            $GTX_TRANSFAC_TFBS, $NCI60, $REPEAT_MASKER, $ENC_HMM, $ENC_TFBS, $SIMPLE_REPEAT, $SEGMENT_DUPS, $ExAC_Freq, $ExAC_AFR, $ExAC_AMR, $ExAC_EAS,
                            $ExAC_FIN, $ExAC_NFE, $ExAC_OTH, $ExAC_SAS, $WELLDERLY_AF, $WELLDERLY_AN,
                            $GME_AF, $GME_NWA, $GME_NEA, $GME_AP, $GME_Israel, $GME_SD, $GME_TP, $GME_CA,
			    $WGSCADD, $WGSCADD_Phred, $WGSDANN, $dbscSNV_ADA_SCORE, $dbscSNV_RF_SCORE, $HRCR1_AF, $HRCR1_AC, $HRCR1_AN, $HRCR1_nonKG_AF, $HRCR1_nonKG_AC, 
                            $HRCR1_nonKG_AN, $WGSFATHMM_coding, $WGSFATHMM_noncoding, $GWAVA_region_score, $GWAVA_tss_score, $ExAC_ALL_nonpsych, $ExAC_ALL_nontcga, 
                            $KAVIAR_AF, $KAVIAR_AC, $KAVIAR_AN, $WGSEIGEN, $RVIS_raw, $RVIS_percent, $ExAC_RVIS_percent, $ExAC_LoF_FDR, $MPC_EXAC, 
                            $UORF_TYPE, $UORF_GENECLASS, $UORF_DISTANCE, $UORF_EFFECT, $UORF_STRENGTH, 
                            $UV_SIGNATURE) {

            ($annot) and $annot =~ s/,$//;
#            ($annot) and $annot =~ s/ /_/g;
#            ($annot) and $annot =~ s/__/_/g;
#            print $annot . "\n";
        }

        # => when number of annotations has to reflect number of alleles, their empty values must be set to 'NIL' for string fields
        foreach my $annot ( $VC, $VT, $SIFT_pred, $Polyphen2_HDIV_pred, $Polyphen2_HVAR_pred, $LRT_pred, $MutationTaster_pred, $MutationAssessor_pred, $FATHMM_pred,
                            $RadialSVM_pred, $LR_pred, $PROVEAN_pred, $FATHMMMKL_coding_pred, $MetaSVM_pred, $MetaLR_pred, $M_CAP_pred, $EIGEN_coding_or_noncoding, $Interpro, $GTEx_V6_gene, $GTEx_V6_tissue, 
                            $CLINVAR, $COSMIC, $ICGC21_ID, $ICGC21_Occurrence, $NCI60, $GTX_GWAS, $GTX_HGMD, $INTERVAR) {
            if (defined $annot) {
                my @values = split /,/, $annot; # / just to avoid a bug in eclipse display...
                if ($annot . "," eq "NIL," x scalar(@values)) {
                    $annot = undef;
#                } else {                   # We have to avoid the case where there is no first value which ends by a leading comma value (e.g. ",bblabla,,blabla")
#                    $annot =~ s/NIL//g;    # -> this is problem because GATK will systematiccaly remove the leading comma, impacting the segregation results !!
		}
            }
        } 

        # Custom ranking of VFT and DA values
        foreach my $annot ($VFT, $DA) {
           next if (!$annot);
            my @values = split /,/, $annot; # / just to avoid a bug in eclipse display...
            foreach my $value (@values) {
               next if ($value !~ /\|/);
                my @unsortedStrings = split /\|/, $value; # / just to avoid a bug in eclipse display...
                my $rA_sortedStrings = customRanking( -input => \@unsortedStrings );
                $value = join "|", @{$rA_sortedStrings};
            }
            my $commaCount = ($annot =~ tr/,//);
            $annot = join ",", @values;
            $annot =~ s/NIL/unknown/g;

        }

        if ($Gene) {
            my @GenesByAllele = split /,/, $Gene; # / just to avoid a bug in eclipse display...
            foreach my $alleleGene (@GenesByAllele) {
                if ($alleleGene =~ /\|/) {
                    my @distinctGenes = split /\|/, $alleleGene; # / just to avoid a bug in eclipse display...
                    foreach my $distinctGene (@distinctGenes) {
                        if ($distinctGene =~ /.+\:.+/) {
                            my @geneTab = split /\:/, $distinctGene; # / just to avoid a bug in eclipse display...
                            $distinctGene = $geneTab[0];
                        }
                    }
                    my %distinctGenes = map { $_ => 1 } sort @distinctGenes;
                    $alleleGene = join "|", sort keys %distinctGenes;
                } else {
                    if ($alleleGene =~ /.+\:.+/) {
                        my @geneTab = split /\:/, $alleleGene; # / just to avoid a bug in eclipse display...
                        $alleleGene = $geneTab[0];
                    }
                }
            }
            $Gene = join ",", sort @GenesByAllele;
        }

        # => when number of annotations has to reflects number of alleles, their empty values must be set to '-1' for float (or interger) fields
        foreach my $annot ( $oneKG, $CG46, $AFR, $AMR, $EUR, $EAS, $SAS, $SIFT_score, $SIFT_rank, $Polyphen2_HDIV_score, $Polyphen2_HDIV_rank, $Polyphen2_HVAR_score, $Polyphen2_HVAR_rank, 
                            $LRT_score, $LRT_rank, $MutationTaster_score, $MutationTaster_rank, $EIGEN_raw, $EIGEN_PC_raw,
                            $MutationAssessor_score, $MutationAssessor_rank, $FATHMM_score, $FATHMM_rank, $RadialSVM_score, $LR_score, $PROVEAN_score, $PROVEAN_rank, 
                            $FATHMMMKL_coding_score, $FATHMMMKL_coding_rank, $MetaSVM_score, $MetaSVM_rank, $MetaLR_score, $MetaLR_rank, $M_CAP_score, $M_CAP_rank,
                            $INTEGRATED_fCons_score, $INTEGRATED_conf_value, $INTEGRATED_fCons_rank, $REVEL,
                            $phyloP100way_vertebrate_score, $phyloP100way_vertebrate_rank, $phyloP20way_mammalian_score, $phyloP20way_mammalian_rank, 
                            $phastCons100way_vertebrate_score, $phastCons100way_vertebrate_rank, $phastCons20way_mammalian_score, $phastCons20way_mammalian_rank,
                            $VEST3_score, $VEST3_rank, $CADD_raw, $CADD_phred, $CADD_rank, $GERP_GT2_WG, $GERP_RS_dbnsfp, $GERP_rank_dbnsfp, $ExAC_Freq, $DANN_score, $DANN_rank, 
                            $GenoCanyon_score, $GenoCanyon_rank, $GME_AF, $GME_NWA, $GME_NEA, $GME_AP, $GME_Israel, $GME_SD, $GME_TP, $GME_CA, 
                            $ExAC_AFR, $ExAC_AMR, $ExAC_EAS, $ExAC_FIN, $ExAC_NFE, $ExAC_OTH, $ExAC_SAS, $SiPhy_29way_logOdds_score, $SiPhy_29way_logOdds_rank, $ESP6500, $ESP6500AA, $ESP6500EA,
                            $gnomAD_exome_ALL, $gnomAD_exome_AFR, $gnomAD_exome_AMR, $gnomAD_exome_ASJ, $gnomAD_exome_EAS, $gnomAD_exome_FIN, $gnomAD_exome_NFE, $gnomAD_exome_OTH, $gnomAD_exome_SAS,
                            $gnomAD_genome_ALL, $gnomAD_genome_AFR, $gnomAD_genome_AMR, $gnomAD_genome_ASJ, $gnomAD_genome_EAS, $gnomAD_genome_FIN, $gnomAD_genome_NFE, $gnomAD_genome_OTH, 
                            $GoNL, $WELLDERLY_AF, $WELLDERLY_AN, $WGSCADD, $WGSCADD_Phred, $WGSDANN, $dbscSNV_ADA_SCORE , $dbscSNV_RF_SCORE  , $HRCR1_AF, $HRCR1_AC, 
                            $HRCR1_AN, $HRCR1_nonKG_AF, $HRCR1_nonKG_AC, $HRCR1_nonKG_AN, $WGSFATHMM_coding, $WGSFATHMM_noncoding, $GWAVA_region_score, 
                            $GWAVA_tss_score, $ExAC_ALL_nonpsych, $ExAC_ALL_nontcga, $KAVIAR_AF, $KAVIAR_AC, $KAVIAR_AN, $WGSEIGEN, $EBI3222PAK, $MPC_EXAC, 
                            $UORF_TYPE, $UORF_GENECLASS, $UORF_DISTANCE, $UORF_EFFECT, $UORF_STRENGTH, 
                            $UV_SIGNATURE){
            if (defined $annot && $annot ne "") {
#print $annot . "\t";
                my @values = split /,/, $annot; # / just to avoid a bug in eclipse display...
#my $tron=eval join '+', @values ;
#my $trin2 = -(scalar(@values));

#print "\"$tron\" \"$trin2\"\n";

                if (-(scalar(@values)) eq eval join '+', @values) {
                    $annot = undef;
#                } else {                   # We have to avoid the case where there is no first value which ends by a leading comma value (e.g. ",bblabla,,blabla")
#                    $annot =~ s/-1//g;     # -> this is problem because GATK will systematiccaly remove the leading comma, impacting the segregation results !!
                }
            }
        }

        my %infoStringHash = (
            'VC'                      => $VC || undef,
            'VFT'                     => $VFT || undef,
            'VT'                      => $VT || undef,
            'Gene'                    => $Gene || undef,
            'DA'                      => $DA || undef,
            'dbSNP_132'               => $dbSNP_132 || undef,
            'dbSNP_138'               => $dbSNP_138 || undef,
            'dbSNP_138Common'         => $dbSNP_138Common || undef,
            'dbSNP_138CommonAlCount'  => $dbSNP_138CommonAlCount || undef,
            'dbSNP_138CommonAlFreq'   => $dbSNP_138CommonAlFreq || undef,
            'dbSNP_138Flagged'        => $dbSNP_138Flagged || undef,
            'dbSNP_138FlaggedAlCount' => $dbSNP_138FlaggedAlCount || undef,
            'dbSNP_138FlaggedAlFreq'  => $dbSNP_138FlaggedAlFreq || undef,
            'dbSNP_142'               => $dbSNP_142 || undef,
            'dbSNP_142Common'         => $dbSNP_142Common || undef,
            'dbSNP_142CommonAlCount'  => $dbSNP_142CommonAlCount || undef,
            'dbSNP_142CommonAlFreq'   => $dbSNP_142CommonAlFreq || undef,
            'dbSNP_142Flagged'        => $dbSNP_142Flagged || undef,
            'dbSNP_142FlaggedAlCount' => $dbSNP_142FlaggedAlCount || undef,
            'dbSNP_142FlaggedAlFreq'  => $dbSNP_142FlaggedAlFreq || undef,
            'KG'                      => $oneKG,
            'KG_AFR'                  => $AFR,
            'KG_AMR'                  => $AMR,
            'KG_EUR'                  => $EUR,
            'KG_EAS'                  => $EAS,
            'KG_SAS'                  => $SAS,
            'ExAC_Freq'               => $ExAC_Freq, 
            'ExAC_AFR'                => $ExAC_AFR,
            'ExAC_AMR'                => $ExAC_AMR,
            'ExAC_EAS'                => $ExAC_EAS,
            'ExAC_FIN'                => $ExAC_FIN,
            'ExAC_NFE'                => $ExAC_NFE,
            'ExAC_OTH'                => $ExAC_OTH,
            'ExAC_SAS'                => $ExAC_SAS,
            'GME_AF'                  => $GME_AF,
            'GME_NWA'                 => $GME_NWA,
            'GME_NEA'                 => $GME_NEA,
            'GME_AP'                  => $GME_AP,
            'GME_Israel'              => $GME_Israel,
            'GME_SD'                  => $GME_SD,
            'GME_TP'                  => $GME_TP,
            'GME_CA'                  => $GME_CA,
            'gnomAD_genome_ALL'       => $gnomAD_genome_ALL,
            'gnomAD_genome_AFR'       => $gnomAD_genome_AFR,
            'gnomAD_genome_AMR'       => $gnomAD_genome_AMR,
            'gnomAD_genome_ASJ'       => $gnomAD_genome_ASJ,
            'gnomAD_genome_EAS'       => $gnomAD_genome_EAS,
            'gnomAD_genome_FIN'       => $gnomAD_genome_FIN,
            'gnomAD_genome_NFE'       => $gnomAD_genome_NFE,
            'gnomAD_genome_OTH'       => $gnomAD_genome_OTH,
            'gnomAD_exome_ALL'        => $gnomAD_exome_ALL,
            'gnomAD_exome_AFR'        => $gnomAD_exome_AFR,
            'gnomAD_exome_AMR'        => $gnomAD_exome_AMR,
            'gnomAD_exome_ASJ'        => $gnomAD_exome_ASJ,
            'gnomAD_exome_EAS'        => $gnomAD_exome_EAS,
            'gnomAD_exome_FIN'        => $gnomAD_exome_FIN,
            'gnomAD_exome_NFE'        => $gnomAD_exome_NFE,
            'gnomAD_exome_OTH'        => $gnomAD_exome_OTH,
            'gnomAD_exome_SAS'        => $gnomAD_exome_SAS,
            'CG46'                    => $CG46,
            'ESP6500'                 => $ESP6500,
            'ESP6500AA'               => $ESP6500AA,
            'ESP6500EA'               => $ESP6500EA,
            'COSMIC'                  => $COSMIC || undef,
            'ICGC21_ID'               => $ICGC21_ID || undef,
            'ICGC21_Occurrence'       => $ICGC21_Occurrence || undef,
            'CLINVAR'                 => $CLINVAR || undef,
            'NCI60'                   => $NCI60 || undef,
            'GoNL'                    => $GoNL || undef,
            'WELLDERLY_AF'            => $WELLDERLY_AF || undef,
            'WELLDERLY_AN'            => $WELLDERLY_AN || undef,
            'SIFT_score'              => $SIFT_score,
            'SIFT_rank'               => $SIFT_rank,
            'SIFT_pred'               => $SIFT_pred,
            'Polyphen2_HDIV_score'    => $Polyphen2_HDIV_score,
            'Polyphen2_HDIV_rank'     => $Polyphen2_HDIV_rank,
            'Polyphen2_HDIV_pred'     => $Polyphen2_HDIV_pred,
            'Polyphen2_HVAR_score'    => $Polyphen2_HVAR_score,
            'Polyphen2_HVAR_rank'     => $Polyphen2_HVAR_rank,
            'Polyphen2_HVAR_pred'     => $Polyphen2_HVAR_pred,
            'LRT_score'               => $LRT_score,
            'LRT_rank'                => $LRT_rank,
            'LRT_pred'                => $LRT_pred,
            'MutationTaster_score'    => $MutationTaster_score,
            'MutationTaster_rank'     => $MutationTaster_rank,
            'MutationTaster_pred'     => $MutationTaster_pred,
            'MutationAssessor_score'  => $MutationAssessor_score,
            'MutationAssessor_rank'   => $MutationAssessor_rank,
            'MutationAssessor_pred'   => $MutationAssessor_pred,
            'FATHMM_score'            => $FATHMM_score,
            'FATHMM_rank'             => $FATHMM_rank,
            'FATHMM_pred'             => $FATHMM_pred,
            'PROVEAN_score'           => $PROVEAN_score,
            'PROVEAN_rank'            => $PROVEAN_rank,
            'PROVEAN_pred'            => $PROVEAN_pred,         
            'DANN_score'              => $DANN_score,
            'DANN_rank'               => $DANN_rank, 
            'FATHMMMKL_coding_score'  => $FATHMMMKL_coding_score,
            'FATHMMMKL_coding_rank'   => $FATHMMMKL_coding_rank,
            'FATHMMMKL_coding_pred'   => $FATHMMMKL_coding_pred,
            'INTEGRATED_fCons_score'  => $INTEGRATED_fCons_score,
            'INTEGRATED_fCons_rank'   => $INTEGRATED_fCons_rank,
            'INTEGRATED_conf_value'   => $INTEGRATED_conf_value,
            'phyloP100way_vertebrate_score' => $phyloP100way_vertebrate_score,
            'phyloP100way_vertebrate_rank'  => $phyloP100way_vertebrate_rank,
            'phyloP20way_mammalian_score'   => $phyloP20way_mammalian_score,
            'phyloP20way_mammalian_rank'    => $phyloP20way_mammalian_rank,
            'phastCons100way_vertebrate_score' => $phastCons100way_vertebrate_score,
            'phastCons100way_vertebrate_rank'  => $phastCons100way_vertebrate_rank,
            'phastCons20way_mammalian_score'   => $phastCons20way_mammalian_score,
            'phastCons20way_mammalian_rank'    => $phastCons20way_mammalian_rank,
            'MetaSVM_score'           => $MetaSVM_score,
            'MetaSVM_rank'            => $MetaSVM_rank,
            'MetaSVM_pred'            => $MetaSVM_pred,
            'MetaLR_score'            => $MetaLR_score,
            'MetaLR_rank'             => $MetaLR_rank,
            'MetaLR_pred'             => $MetaLR_pred,
            'VEST3_score'             => $VEST3_score,
            'VEST3_rank'              => $VEST3_rank,
            'M_CAP_score'             => $M_CAP_score,
            'M_CAP_rank'              => $M_CAP_rank,
            'M_CAP_pred'              => $M_CAP_pred,
            'EIGEN_coding_or_noncoding' => $EIGEN_coding_or_noncoding,
            'EIGEN_raw'               => $EIGEN_raw,
            'EIGEN_PC_raw'            => $EIGEN_PC_raw,
            'GenoCanyon_score'        => $GenoCanyon_score,
            'GenoCanyon_rank'         => $GenoCanyon_rank,
            'CADD_raw'                => $CADD_raw,
            'CADD_rank'               => $CADD_rank,
            'CADD_phred'              => $CADD_phred,
            'GERP_RS_dbnsfp'                 => $GERP_RS_dbnsfp,
            'GERP_rank_dbnsfp'               => $GERP_rank_dbnsfp,
            'GERP_GT2_WG'               => $GERP_GT2_WG,
            'GERP_ELEM_WG'              => $GERP_ELEM_WG || undef,
            'SiPhy_29way_logOdds_score'     => $SiPhy_29way_logOdds_score,
            'SiPhy_29way_logOdds_rank'      => $SiPhy_29way_logOdds_rank,
            'Interpro'                => $Interpro,
            'GTEx_V6_gene'            => $GTEx_V6_gene,
            'GTEx_V6_tissue'          => $GTEx_V6_tissue,
            'MIRNA'                   => $MIRNA || undef,
            'MIRNATARGET'             => $MIRNATARGET || undef,
            'RNASEQ_CUFF'             => $RNASEQ_CUFF || undef,
            'RNASEQ_SCRI'             => $RNASEQ_SCRI || undef,
            'GTX_CHIPSEQ'             => $GTX_CHIPSEQ || undef,
            'GTX_CPG'                 => $GTX_CPG || undef,
            'GTX_DISEASE'             => $GTX_DISEASE || undef,
            'GTX_DNASE'               => $GTX_DNASE || undef,
            'GTX_DRUG'                => $GTX_DRUG || undef,
            'GTX_GWAS'                => $GTX_GWAS || undef,
            'GTX_HGMD'                => $GTX_HGMD || undef,
            'GTX_HGMD_DISGENE'        => $GTX_HGMD_DISGENE || undef,
            'GTX_HGMDIMPUTED'         => $GTX_HGMDIMPUTED || undef,
            'GTX_MICROSAT'            => $GTX_MICROSAT || undef,
            'GTX_MIRNA'               => $GTX_MIRNA || undef,
            'GTX_OMIM'                => $GTX_OMIM || undef,
            'GTX_ORPHA'               => $GTX_ORPHA || undef,
            'GTX_PATH'                => $GTX_PATH || undef,
            'GTX_PGMD'                => $GTX_PGMD || undef,
            'GTX_PTMS'                => $GTX_PTMS || undef,
            'GTX_TRANSFAC_TFBS'       => $GTX_TRANSFAC_TFBS || undef,
            'GTX_TSS'                 => $GTX_TSS || undef,
            'GTX_VISTA'               => $GTX_VISTA || undef,
            'ENC_HMM'                 => $ENC_HMM || undef,
            'ENC_TFBS'                => $ENC_TFBS || undef,
            'SIMPLE_REPEAT'           => $SIMPLE_REPEAT || undef,
            'SEGMENT_DUPS'            => $SEGMENT_DUPS || undef,
            'REPEAT_MASKER'           => $REPEAT_MASKER || undef,
	    'WGSCADD'                 => $WGSCADD, 
            'WGSCADD_Phred'           => $WGSCADD_Phred, 
            'WGSDANN'                 => $WGSDANN, 
            'dbscSNV_ADA_SCORE'       => $dbscSNV_ADA_SCORE,
            'dbscSNV_RF_SCORE'        => $dbscSNV_RF_SCORE, 
            'HRCR1_AF'                => $HRCR1_AF, 
	    'HRCR1_AC'                => $HRCR1_AC,
            'HRCR1_AN'                => $HRCR1_AN,
            'HRCR1_nonKG_AF'          => $HRCR1_nonKG_AF,
            'HRCR1_nonKG_AC'          => $HRCR1_nonKG_AC,
            'HRCR1_nonKG_AN'          => $HRCR1_nonKG_AN,
            'WGSFATHMM_coding'        => $WGSFATHMM_coding, 
            'WGSFATHMM_noncoding'     => $WGSFATHMM_noncoding,
            'GWAVA_region_score'      => $GWAVA_region_score,
            'GWAVA_tss_score'         => $GWAVA_tss_score,
            'ExAC_ALL_nonpsych'       => $ExAC_ALL_nonpsych,
	    'ExAC_ALL_nontcga'        => $ExAC_ALL_nontcga,
            'KAVIAR_AF'               => $KAVIAR_AF,            
            'KAVIAR_AC'               => $KAVIAR_AC,
            'KAVIAR_AN'               => $KAVIAR_AN,
            'WGSEIGEN'                => $WGSEIGEN,
            'RVIS_raw'                => $RVIS_raw || undef,
            'RVIS_percent'            => $RVIS_percent || undef,
            'ExAC_RVIS_percent'       => $ExAC_RVIS_percent || undef,
            'ExAC_LoF_FDR'            => $ExAC_LoF_FDR || undef,
            'REVEL'                   => $REVEL,
            'INTERVAR'                => $INTERVAR,
            'EBI3222PAK'              => $EBI3222PAK,
            'MPC_EXAC'                => $MPC_EXAC,
            'uORF_type'               => $UORF_TYPE,
            'uORF_geneClass'          => $UORF_GENECLASS,
            'uORF_distance'           => $UORF_DISTANCE,
            'uORF_effect'             => $UORF_EFFECT,
            'uORF_kozakStrength'      => $UORF_STRENGTH,
            'UV_signature'            => $UV_SIGNATURE,
        );

        # Manage info fields which have no value associated <=> just a boolean flag
        ($infoStringHash{dbSNP_138Common})  and $infoStringHash{dbSNP_138Common}  = "";
        ($infoStringHash{dbSNP_138Flagged}) and $infoStringHash{dbSNP_138Flagged} = "";
        ($infoStringHash{dbSNP_142Common})  and $infoStringHash{dbSNP_142Common}  = "";
        ($infoStringHash{dbSNP_142Flagged}) and $infoStringHash{dbSNP_142Flagged} = "";

        foreach my $infoTag (qw( gene snp132 snp138 snp138common snp138flagged snp142 snp142common snp142flagged 1kg 1kg_afr 1kg_amr 1kg_eur 1kg_eas 1kg_sas esp6500 esp6500aa esp6500ea nci60 cosmic icgc21
                                 exac exac_afr exac_amr exac_eas exac_fin exac_nfe exac_oth exac_sas cg gme gnomad_exome gnomad_genome clinvar gonl mirna mirnatarget rnaseq_cuff rnaseq_scri rvis_exac
                                 dbnsfp sift pp2 lrt mt ma fathmm rsvm lr vest cadd gerp gerp++gt2 gerp++elem phylop siphy mcap eigen genocanyon interpro gtex enc_hmm enc_tfbs simple_repeat
                                 segment_dups rmsk wellderly gtx_chipseq gtx_cpg gtx_disease gtx_dnase gtx_drug gtx_gwas gtx_hgmd gtx_hgmd_disgene gtx_hgmdimputed gtx_microsat gtx_mirna
                                 gtx_omim gtx_orpha gtx_path gtx_pgmd gtx_ptms gtx_transfac gtx_tss gtx_vista ebi3222pak mpc_exac uorf uv_signature)) {

            if (!doIt(-elem => $infoTag)) {
                if ($infoTag eq "gene") {
                    delete $infoStringHash{Gene};
                    delete $infoStringHash{VFT};
                    delete $infoStringHash{VT};
                    delete $infoStringHash{VC};
                    delete $infoStringHash{DA};
                } elsif ($infoTag eq "snp132") {
                    delete $infoStringHash{dbSNP_132};
                } elsif ($infoTag eq "snp138") {
                    delete $infoStringHash{dbSNP_138};
                } elsif ($infoTag eq "snp138common") {
                    delete $infoStringHash{dbSNP_138Common};
                    delete $infoStringHash{dbSNP_138CommonAlCount};
                    delete $infoStringHash{dbSNP_138CommonAlFreq};
                }  elsif ($infoTag eq "snp138flagged") {
                    delete $infoStringHash{dbSNP_138Flagged};
                    delete $infoStringHash{dbSNP_138FlaggedAlCount};
                    delete $infoStringHash{dbSNP_138FlaggedAlFreq};
                } elsif ($infoTag eq "snp142") {
                    delete $infoStringHash{dbSNP_142};
                } elsif ($infoTag eq "snp142common") {
                    delete $infoStringHash{dbSNP_142Common};
                    delete $infoStringHash{dbSNP_142CommonAlCount};
                    delete $infoStringHash{dbSNP_142CommonAlFreq};
                }  elsif ($infoTag eq "snp142flagged") {
                    delete $infoStringHash{dbSNP_142Flagged};
                    delete $infoStringHash{dbSNP_142FlaggedAlCount};
                    delete $infoStringHash{dbSNP_142FlaggedAlFreq};
                } elsif ($infoTag eq "1kg_eas") {
                    delete $infoStringHash{KG_EAS};
                } elsif ($infoTag eq "1kg_sas") {
                    delete $infoStringHash{KG_SAS};
                } elsif ($infoTag eq "1kg_afr") {
                    delete $infoStringHash{KG_AFR};
                } elsif ($infoTag eq "1kg_amr") {
                    delete $infoStringHash{KG_AMR};
                } elsif ($infoTag eq "1kg_eur") {
                    delete $infoStringHash{KG_EUR};
                } elsif ($infoTag eq "gonl") {
                    delete $infoStringHash{GoNL};
                } elsif ($infoTag eq "wellderly") {
                    delete $infoStringHash{WELLDERLY_AF};
                    delete $infoStringHash{WELLDERLY_AN};
                } elsif ($infoTag eq "exac") {
                    delete $infoStringHash{ExAC_Freq};
                    delete $infoStringHash{ExAC_AFR};
                    delete $infoStringHash{ExAC_AMR};
                    delete $infoStringHash{ExAC_EAS};
                    delete $infoStringHash{ExAC_FIN};
                    delete $infoStringHash{ExAC_NFE};
                    delete $infoStringHash{ExAC_OTH};
                    delete $infoStringHash{ExAC_SAS};
                } elsif ($infoTag eq "exac_afr") {
                    delete $infoStringHash{ExAC_AFR};
                } elsif ($infoTag eq "exac_amr") {
                    delete $infoStringHash{ExAC_AMR};
                } elsif ($infoTag eq "exac_eas") {
                    delete $infoStringHash{ExAC_EAS};
                } elsif ($infoTag eq "exac_fin") {
                    delete $infoStringHash{ExAC_FIN};
                } elsif ($infoTag eq "exac_nfe") {
                    delete $infoStringHash{ExAC_NFE};
                } elsif ($infoTag eq "exac_oth") {
                    delete $infoStringHash{ExAC_OTH};
                } elsif ($infoTag eq "exac_sas") {
                    delete $infoStringHash{ExAC_SAS};
                } elsif ($infoTag eq "gme") {
                    delete $infoStringHash{"GME_AF"};
                    delete $infoStringHash{"GME_NWA"};
                    delete $infoStringHash{"GME_NEA"};
                    delete $infoStringHash{"GME_AP"};
                    delete $infoStringHash{"GME_Israel"};
                    delete $infoStringHash{"GME_SD"};
                    delete $infoStringHash{"GME_TP"};
                    delete $infoStringHash{"GME_CA"};
                } elsif ($infoTag eq "gnomad_exome") {
                    delete $infoStringHash{"gnomAD_exome_ALL"};
                    delete $infoStringHash{"gnomAD_exome_AFR"};
                    delete $infoStringHash{"gnomAD_exome_AMR"};
                    delete $infoStringHash{"gnomAD_exome_ASJ"};
                    delete $infoStringHash{"gnomAD_exome_EAS"};
                    delete $infoStringHash{"gnomAD_exome_FIN"};
                    delete $infoStringHash{"gnomAD_exome_NFE"};
                    delete $infoStringHash{"gnomAD_exome_OTH"};
                    delete $infoStringHash{"gnomAD_exome_SAS"};
                } elsif ($infoTag eq "gnomad_genome") {
                    delete $infoStringHash{"gnomAD_genome_ALL"};
                    delete $infoStringHash{"gnomAD_genome_AFR"};
                    delete $infoStringHash{"gnomAD_genome_AMR"};
                    delete $infoStringHash{"gnomAD_genome_ASJ"};
                    delete $infoStringHash{"gnomAD_genome_EAS"};
                    delete $infoStringHash{"gnomAD_genome_FIN"};
                    delete $infoStringHash{"gnomAD_genome_NFE"};
                    delete $infoStringHash{"gnomAD_genome_OTH"};
                } elsif ($infoTag eq "dbnsfp") {
		    (!doIt(-elem => "sift"))     and delete $infoStringHash{SIFT_score};
                    (!doIt(-elem => "sift"))     and delete $infoStringHash{SIFT_pred};
                    (!doIt(-elem => "sift"))     and delete $infoStringHash{SIFT_rank};
                    (!doIt(-elem => "pp2"))      and delete $infoStringHash{Polyphen2_HDIV_score};
                    (!doIt(-elem => "pp2"))      and delete $infoStringHash{Polyphen2_HDIV_rank};
                    (!doIt(-elem => "pp2"))      and delete $infoStringHash{Polyphen2_HDIV_pred};
                    (!doIt(-elem => "pp2"))      and delete $infoStringHash{Polyphen2_HVAR_score};
                    (!doIt(-elem => "pp2"))      and delete $infoStringHash{Polyphen2_HVAR_rank};
                    (!doIt(-elem => "pp2"))      and delete $infoStringHash{Polyphen2_HVAR_pred};
                    (!doIt(-elem => "lrt"))      and delete $infoStringHash{LRT_score};
                    (!doIt(-elem => "lrt"))      and delete $infoStringHash{LRT_rank};
                    (!doIt(-elem => "lrt"))      and delete $infoStringHash{LRT_pred};
                    (!doIt(-elem => "mt"))       and delete $infoStringHash{MutationTaster_score};
                    (!doIt(-elem => "mt"))       and delete $infoStringHash{MutationTaster_rank};
                    (!doIt(-elem => "mt"))       and delete $infoStringHash{MutationTaster_pred};
                    (!doIt(-elem => "ma"))       and delete $infoStringHash{MutationAssessor_score};
                    (!doIt(-elem => "ma"))       and delete $infoStringHash{MutationAssessor_rank};
                    (!doIt(-elem => "ma"))       and delete $infoStringHash{MutationAssessor_pred};
                    (!doIt(-elem => "fathmm"))   and delete $infoStringHash{FATHMM_score};
                    (!doIt(-elem => "fathmm"))   and delete $infoStringHash{FATHMM_rank};
                    (!doIt(-elem => "fathmm"))   and delete $infoStringHash{FATHMM_pred};
		    (!doIt(-elem => "provean"))  and delete $infoStringHash{PROVEAN_score};
                    (!doIt(-elem => "provean"))  and delete $infoStringHash{PROVEAN_rank};
		    (!doIt(-elem => "provean"))  and delete $infoStringHash{PROVEAN_pred};
                    (!doIt(-elem => "vest"))     and delete $infoStringHash{VEST3_score};
                    (!doIt(-elem => "vest"))     and delete $infoStringHash{VEST3_rank};
                    (!doIt(-elem => "mcap"))     and delete $infoStringHash{M_CAP_score};
                    (!doIt(-elem => "mcap"))     and delete $infoStringHash{M_CAP_rank};
                    (!doIt(-elem => "mcap"))     and delete $infoStringHash{M_CAP_pred};
                    (!doIt(-elem => "cadd"))     and delete $infoStringHash{CADD_raw};
                    (!doIt(-elem => "cadd"))     and delete $infoStringHash{CADD_rank};
                    (!doIt(-elem => "cadd"))     and delete $infoStringHash{CADD_phred};
		    (!doIt(-elem => "dann"))     and delete $infoStringHash{DANN_score};
                    (!doIt(-elem => "dann"))     and delete $infoStringHash{DANN_rank};
		    (!doIt(-elem => "fathmmkl")) and delete $infoStringHash{FATHMMMKL_coding_score};
                    (!doIt(-elem => "fathmmkl")) and delete $infoStringHash{FATHMMMKL_coding_rank};
		    (!doIt(-elem => "fathmmkl")) and delete $infoStringHash{FATHMMMKL_coding_pred};
                    (!doIt(-elem => "eigen"))    and delete $infoStringHash{EIGEN_coding_or_noncoding};
                    (!doIt(-elem => "eigen"))    and delete $infoStringHash{EIGEN_raw};
                    (!doIt(-elem => "eigen"))    and delete $infoStringHash{EIGEN_PC_raw};
                    (!doIt(-elem => "genocanyon")) and delete $infoStringHash{GenoCanyon_score};
                    (!doIt(-elem => "genocanyon")) and delete $infoStringHash{GenoCanyon_rank};
                    (!doIt(-elem => "msvm"))     and delete $infoStringHash{MetaSVM_score};
                    (!doIt(-elem => "msvm"))     and delete $infoStringHash{MetaSVM_rank};
                    (!doIt(-elem => "msvm"))     and delete $infoStringHash{MetaSVM_pred};
                    (!doIt(-elem => "mlr"))      and delete $infoStringHash{MetaLR_score};
                    (!doIt(-elem => "mlr"))      and delete $infoStringHash{MetaLR_rank};
                    (!doIt(-elem => "mlr"))      and delete $infoStringHash{MetaLR_pred};
		    (!doIt(-elem => "intgr"))    and delete $infoStringHash{INTEGRATED_fCons_score};
                    (!doIt(-elem => "intgr"))    and delete $infoStringHash{INTEGRATED_fCons_rank};
		    (!doIt(-elem => "intgr"))    and delete $infoStringHash{INTEGRATED_conf_value};
                    (!doIt(-elem => "gerp"))     and delete $infoStringHash{GERP_RS_dbnsfp};
                    (!doIt(-elem => "gerp"))     and delete $infoStringHash{GERP_rank_dbnsfp};
                    (!doIt(-elem => "phylop"))   and delete $infoStringHash{phyloP100way_vertebrate_score};
                    (!doIt(-elem => "phylop"))   and delete $infoStringHash{phyloP100way_vertebrate_rank};
                    (!doIt(-elem => "phylop"))   and delete $infoStringHash{phyloP20way_mammalian_score};
                    (!doIt(-elem => "phylop"))   and delete $infoStringHash{phyloP20way_mammalian_rank};
		    (!doIt(-elem => "phcons"))   and delete $infoStringHash{phastCons100way_vertebrate_score};
		    (!doIt(-elem => "phcons"))   and delete $infoStringHash{phastCons100way_vertebrate_rank};
                    (!doIt(-elem => "phcons"))   and delete $infoStringHash{phastCons20way_mammalian_score};
                    (!doIt(-elem => "phcons"))   and delete $infoStringHash{phastCons20way_mammalian_rank};
                    (!doIt(-elem => "siphy"))    and delete $infoStringHash{SiPhy_29way_logOdds_score};
                    (!doIt(-elem => "siphy"))    and delete $infoStringHash{SiPhy_29way_logOdds_rank};
                    (!doIt(-elem => "interpro")) and delete $infoStringHash{Interpro};
                    (!doIt(-elem => "gtex"))     and delete $infoStringHash{GTEx_V6_gene};
                    (!doIt(-elem => "gtex"))     and delete $infoStringHash{GTEx_V6_tissue};
                } elsif ($infoTag eq "sift") {
                    delete $infoStringHash{SIFT_score};
                    delete $infoStringHash{SIFT_rank};
                    delete $infoStringHash{SIFT_pred};
                } elsif ($infoTag eq "pp2") {
                    delete $infoStringHash{Polyphen2_HDIV_score};
                    delete $infoStringHash{Polyphen2_HDIV_rank};
                    delete $infoStringHash{Polyphen2_HDIV_pred};
                    delete $infoStringHash{Polyphen2_HVAR_score};
                    delete $infoStringHash{Polyphen2_HVAR_rank};
                    delete $infoStringHash{Polyphen2_HVAR_pred};
                } elsif ($infoTag eq "lrt") {
                    delete $infoStringHash{LRT_score};
                    delete $infoStringHash{LRT_rank};
                    delete $infoStringHash{LRT_pred};
                } elsif ($infoTag eq "mt") {
                    delete $infoStringHash{MutationTaster_score};
                    delete $infoStringHash{MutationTaster_rank};
                    delete $infoStringHash{MutationTaster_pred};
                } elsif ($infoTag eq "ma") {
                    delete $infoStringHash{MutationAssessor_score};
                    delete $infoStringHash{MutationAssessor_rank};
                    delete $infoStringHash{MutationAssessor_pred};
                } elsif ($infoTag eq "fathmm") {
                    delete $infoStringHash{FATHMM_score};
                    delete $infoStringHash{FATHMM_rank};
                    delete $infoStringHash{FATHMM_pred};
		} elsif ($infoTag eq "provean") {
		    delete $infoStringHash{PROVEAN_score};
                    delete $infoStringHash{PROVEAN_rank};
		    delete $infoStringHash{PROVEAN_pred};
                } elsif ($infoTag eq "vest") {
                    delete $infoStringHash{VEST3_score};
                    delete $infoStringHash{VEST3_rank};
                } elsif ($infoTag eq "mcap") {
                    delete $infoStringHash{M_CAP_score};
                    delete $infoStringHash{M_CAP_rank};
                    delete $infoStringHash{M_CAP_pred};
                } elsif ($infoTag eq "cadd") {
                    delete $infoStringHash{CADD_raw};
                    delete $infoStringHash{CADD_rank};
                    delete $infoStringHash{CADD_phred};
                } elsif ($infoTag eq "dann") {
		    delete $infoStringHash{DANN_score};
                    delete $infoStringHash{DANN_rank};
		} elsif ($infoTag eq "fathmmkl") {
		    delete $infoStringHash{FATHMMMKL_coding_score};
                    delete $infoStringHash{FATHMMMKL_coding_rank};
		    delete $infoStringHash{FATHMMMKL_coding_pred};
                } elsif ($infoTag eq "eigen") {
                    delete $infoStringHash{EIGEN_coding_or_noncoding};
                    delete $infoStringHash{EIGEN_raw};
                    delete $infoStringHash{EIGEN_PC_raw};
                } elsif ($infoTag eq "genocanyon") {
                    delete $infoStringHash{GenoCanyon_score};
                    delete $infoStringHash{GenoCanyon_rank};
		} elsif ($infoTag eq "msvm") {
                    delete $infoStringHash{MetaSVM_score};
                    delete $infoStringHash{MetaSVM_rank};
                    delete $infoStringHash{MetaSVM_pred};
                } elsif ($infoTag eq "mlr") {
                    delete $infoStringHash{MetaLR_score};
                    delete $infoStringHash{MetaLR_rank};
                    delete $infoStringHash{MetaLR_pred};
		} elsif ($infoTag eq "intgr") {
		    delete $infoStringHash{INTEGRATED_fCons_score};
                    delete $infoStringHash{INTEGRATED_fCons_rank};
		    delete $infoStringHash{INTEGRATED_conf_value};
                } elsif ($infoTag eq "gerp") {
                    delete $infoStringHash{GERP_RS_dbnsfp};
                    delete $infoStringHash{GERP_rank_dbnsfp};
                } elsif ($infoTag eq "phylop") {
                    delete $infoStringHash{phyloP100way_vertebrate_score};
                    delete $infoStringHash{phyloP100way_vertebrate_rank};
                    delete $infoStringHash{phyloP20way_mammalian_score};
                    delete $infoStringHash{phyloP20way_mammalian_rank};
                } elsif ($infoTag eq "siphy") {
                    delete $infoStringHash{SiPhy_29way_logOdds_score};
                    delete $infoStringHash{SiPhy_29way_logOdds_rank};
		} elsif ($infoTag eq "phcons") {
		     delete $infoStringHash{phastCons100way_vertebrate_score};
                     delete $infoStringHash{phastCons100way_vertebrate_rank};
                     delete $infoStringHash{phastCons20way_mammalian_score};
                     delete $infoStringHash{phastCons20way_mammalian_rank};
                } elsif ($infoTag eq "interpro") {
                     delete $infoStringHash{Interpro};
                } elsif ($infoTag eq "gtex") {
                     delete $infoStringHash{GTEx_V6_gene};
                     delete $infoStringHash{GTEx_V6_tissue};
                } elsif ($infoTag eq "rvis_exac") {
                     delete $infoStringHash{RVIS_raw};
                     delete $infoStringHash{RVIS_percent};
                     delete $infoStringHash{ExAC_RVIS_percent};
                     delete $infoStringHash{ExAC_LoF_FDR};
                } elsif ($infoTag eq "ebi3222pak") {
                     delete $infoStringHash{EBI3222PAK};
                } elsif ($infoTag eq "mpc_exac") {
                     delete $infoStringHash{MPC_EXAC};
                } elsif ($infoTag eq "uorf") {
                     delete $infoStringHash{uORF_type};
                     delete $infoStringHash{uORF_geneClass};
                     delete $infoStringHash{uORF_distance};
                     delete $infoStringHash{uORF_effect};
                     delete $infoStringHash{uORF_kozakStrength};
                } elsif ($infoTag eq "uv_signature") {
                     delete $infoStringHash{UV_signature};
                } else {
                    delete $infoStringHash{uc($infoTag)};
                }
            }
        }
        $rA_x->[7] = $rO_vcf->add_info_field($rA_x->[7], %infoStringHash);
        $rA_x->[7] = join ";", sort split /;/, $rA_x->[7];

        print RESULTS $rO_vcf->format_line($rA_x);
    }
    $rO_vcf->close();
    close RESULTS;
}

sub doIt {
    my %arg = @_;
    my $element = $arg{-elem};
    (!$EXCLUDE && !$INCLUDE) and return 1;
    if ($INCLUDE) {
        if ($INCLUDE eq "all") {
            (inExcTab(-elem => $element)) ? return 0 : return 1;
        } elsif (inIncTab(-elem => $element)) {
             return 1;
        } else {
            return 0;
        }
    }
    if ($EXCLUDE) {
        if ($EXCLUDE eq "all") {
            ($element eq "gene") and return 1;
            (inIncTab(-elem => $element)) ? return 1 : return 0;
        } elsif (inExcTab(-elem => $element)) {
             return 0;
        } else {
            return 1;
        }
    }
}

sub customRanking {
    my %arg = @_;
    my $rA_input = $arg{-input};

    my @input = @{$rA_input};

    my @FUNCTION_RANKING = qw(exonic exonic_splicing intronic_splicing intronic ncRNA_exonic ncRNA_splicing ncRNA_intronic UTR3_splicing UTR5_splicing UTR3 UTR5 downstream upstream intergenic);
    my %FUNCTION_RANKING_MAP = map { $FUNCTION_RANKING[$_] => $_ } 0 .. $#FUNCTION_RANKING;
    my $RANKING_PATTERN = join '|', @FUNCTION_RANKING;

    my @sorted = reverse sort {
        my ($x, $y) = map /($RANKING_PATTERN)/, $a, $b;
        $FUNCTION_RANKING_MAP{$y} <=> $FUNCTION_RANKING_MAP{$x};
    } @input;
    return \@sorted;
}

sub inIncTab {
    my %arg = @_;
    my $element = $arg{-elem};

    my %hash = map { lc($_) => 1 } @INCLUDETAB;

    return $hash{$element};
}

sub inExcTab {
    my %arg = @_;
    my $element = $arg{-elem};

    my %hash = map { lc($_) => 1 } @EXCLUDETAB;

    return $hash{$element};
}

sub GetCurrentTime {
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    $month++;
    return sprintf "%04d-%02d-%02d_%02dh%02dm%02ds", ($year, $month, $dayOfMonth, $hour, $minute, $second);
}

sub printHeader {
    my %arg = @_;
    my $relultFH = $arg{-fh};

    open(HEAD, "<$HEADER") || die "Could not open file $HEADER...\n$!\n";

    while (my $line = <HEAD>) {
        chomp $line;
        print $relultFH $line, "\n";
    }
}

sub printerr {
    print STDERR @_;
    print LOG @_;
}

=head1 SYNOPSIS

 vcf2Annovar.pl [arguments]

 Mandatory arguments to control input and output files
        -i,  --infile     <file_path>    path of the input file(s)
        -b,  --buildver   <string>       genome build version (hg18 (or 18), hg19 (or 19), b37 or v37)
        -vc, --vcaller    <string>       name of the software (algorithm) that did generate the input files (varscan, diBayes, smallIndel(s), mpileup, gatk, dindel, custom...)
        -id, --identifier <string>       specify the identificator (e.g. sample) on which the script will run (used to create proper working files...)

 Optional arguments:
        -o,   --outfile   <file_path>    path of the output file(s) (default: current user directory, filename generated based on variant_caller, identifier...)
        -inc, --include   <string>       database filter(s) to include in analysis (e.g snp or snp,1kg,cg), cannot be combined with --exclude unless --exclude is set to 'all'
        -e,   --exclude   <string>       database filter(s) to exclude from analysis (e.g snp or snp,1kg,cg), cannot be combined with --include unless --include is set to 'all'
        -rf,  --rfolder   <string>       run the script in resume mode : path of the folder containing all the previouisly computed files (annovar input file, .dropped files...)
        -h,   --help                     print help message
        -m,   --man                      print complete documentation
        -v,   --verbose                  use verbose output

 Function: annotate a vcf input with Annovar, then outputting result in a new vcf file containing all the requested annotations in the INFO field and an updated vcf header.

 Version: $LastChangedDate: 2011-11-16 15:35:31 -0500 (Wed, 16 Nov 2011) $

=head1 OPTIONS

=over 12

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--infile>

path of the input file (generated by varscan or dibayes or smallindels or etc...) containing all the variants information.
no default value : it has to be provided !

=item B<--outfile>

path of the output results file. By default, the same path as the provided -infile will be used.
each other output (temporary, dropped, filtered, log, etc...) files will be written in the same path as the output file.

=item B<--identifier>

identifier of the current analysis.
e.g. : name of the sample we want to analyse (if so, must correspond to S2D database sample).

=item B<--buildver>

genome build version to use. The build version will be used by ANNOVAR to identify
corresponding database files automatically, for example, when gene-based annotation
is used for hg18 build, ANNOVAR will search for the hg18_refGene.txt file, but if
the hg19 is used as --buildver, ANNOVAR will examine hg19_refGene.txt instead.

=item B<--vcaller>

specifies wich variant caller has been used to genearate the input file.
until now, seven variant callers are supported : varscan, diBayes, smallIndel(s), mpileup, gatk, dindel, custom

=item B<--include>

specifies the database filters which the user wants to include in his analysis, cannot be used with --exclude unless --exclude is set to 'all'.
supported choices are : gene (include gene annotation), snp132 (include dbSNP_132 filter), snp138 (include dbSNP_138 filter), snp138common (include common SNPs from dbSNP_138 filter), snp138flagged (include flagged SNPs from dbSNP_138 filter), snp142 (include dbSNP_142 filter), snp142common (include common SNPs from dbSNP_142 filter), snp142flagged (include flagged SNPs from dbSNP_142 filter), 1kg (include 1000 Genome filter), 1kg_afr (include filter of allele frequency of African population of 1000genomes project), 1kg_amr (include filter of allele frequency of Ad Mixed American population of 1000genomes project), 1kg_eur (include filter of allele frequency of European population of 1000genomes project), 1kg_eas (include filter of allele frequency of East Asian population of 1000genomes project), 1kg_sas (include filter of allele frequency of South Asian population of 1000genomes project), exac (include filter), exac_afr (include filter), exac_amr (include filter), exac_eas (include filter), exac_fin (include filter), exac_nfe (include filter), exac_oth (include filter), exac_sas (include filter), cg (include Complete Genomics filter), gme (include Great Middle East frequencies), dbnsfp (include dbnsfp33a database filter : including whole-exome SIFT scores, PolyPhen2 HDIV scores, PolyPhen2 HVAR scores, LRT scores, MutationTaster scores, MutationAssessor score, FATHMM scores, PROVEAN scores, VEST scores, CADD scores, DANN scores, M-CAP scores, Eigen scores, GenoCanyon scores, Interpro domains, GTex annotations, FATHMMMKL scores, MetaSVM scores, MetaLR scores, INTEGRATED scores, GERP++ scores, PhyloP scores, Phast scores and SiPhy scores from dbnsfp (dbnsfp version 3.3a - 20170221), sift (include SIFT filter), pp2 (include PolyPhen.v2 HumanDiv & HumanVar filters), mt (include Mutation Taster filter), lrt (include LRT filter), phylop (include PhyloP filter), gerp_gt2 (whole-genome gerp++ scores greater than 2), gerp_elem (conserved genomic regions by gerp++), ma (include filter), fathmm (include filter), rsvm (include filter), lr (include filter), vest (include filter), cadd (include filter), gerp (include filter), siphy (include filter), esp6500 (include esp6500 filter), esp6500aa (include esp6500aa filter), esp6500ea (include esp6500ea filter), mirna (include mirna region annotation), mirnatarget (include mirnatarget region annotation), rnaseq_cuff (include Broad Cufflinks RNASeq alignment region annotation), rnaseq_scri (include Broad Scripture RNASeq alignment region annotation), clinvar (include NCBI\'s ClinVar annotation), cosmic (include Sanger's COSMIC (Catalogue of Somatic Mutations in Cancer) somatic disease mutations description region annotation), icgc21 (International Cancer Genome Consortium version 21), nci60 (include NCI-60 allele frequency filter), gtx_chipseq (include GenomeTrax Predicted ChIP-Seq TFBS description region annotation), gtx_cpg (include GenomeTrax CpG Islands description region annotation), gtx_disease (include GenomeTrax Disease associations description region annotation), gtx_dnase (Include GenomeTrax Dnase biding site), gtx_drug (include GenomeTrax Drug Targets description region annotation), gtx_gwas (include GenomeTrax GWAS Catalogue description region annotation), gtx_hgmd (include GenomeTrax HGMD inherited (germ-line) disease mutations description region annotation), gtx_hgmd_disgene (include GenomeTrax HGMD disease genes description region annotation), gtx_hgmdimputed (include GenomeTrax HGMD IMPUTED mutation), gtx_microsat (include GenomeTrax Microsatellite repeats description region annotation), gtx_mirna (include GenomeTrax microRNA sequences description region annotation), gtx_omim (include GenomeTrax NCBI\'s OMIM description region annotation), gtx_orpha (include GenomeTrax rare diseases information), gtx_path (include GenomeTrax Pathways membership description region annotation), gtx_pgmd (include GenomeTrax Pharmacogenomics Variants description region annotation), gtx_ptms (include GenomeTrax Post translational modifications description region annotation), gtx_transfac (include GenomeTrax TRANSFAC experimentally verified TFBS description region annotation), gtx_tss (include GenomeTrax TSSs (transcription start sites) description region annotation), gtx_vista (include vista for tracking predicted damaging effect of a variation occurring in an enhancer site), enc_hmm (include Encode Chromatin States (UCSC track) description region annotation), enc_tfbs (include Encode Transcription Factor Binding Sites (UCSC track) description region annotation), simple_repeat (include Simple Tandem Repeats (possibly imperfect repeats) located by Tandem Repeats Finder (TRF) (UCSC simpleRepeat track) description region annotation), segment_dups (include regions detected as putative genomic duplications (UCSC  genomicSuperDups track) description region annotation), rmsk (include regions detected by RepeatMasker as interspersed repeats and low complexity DNA sequencesinterspersed repeats and low complexity DNA sequences), gonl (include GoNL project frequencies), wgsdann (include dann whole-genome variants score), gwava (include gwava region/tss score), wgscadd (include cadd/cadd_Phred whole genome score), wgsfathmm (include fathmm whole-genome coding/noncoding score), hrcr1 (include 40M variants from 32K samples (HRC_[AF,AC,AN] AND HRC_nonKG_[AF,AC,AN])), dbscsnv11 (include dbscsnv splice site effect prediction), exac_nontcga (include exac nontcga allele frequency), exac_nonpsych (include exac non psych allele frequency), kaviar (include kaviar allele frequency), wgseigen (include wg eigen coding & noncoding socre), rvis_exac (include RVIS (Residual Variation Intolerance Score) score based on ExAC sequencing data), revel (include REVEL score for predicting pathogenicity of missense variants based on a combination of scores from 13 individual tools), intervar (include INTERVAR score for clinical interpretation of missense variants), ebi3222pak (include MAF from EBI dataset for 3222 healthy pakistani controls), mpc_exac (ExAC regional missense constraint score), uorf (5'UTR variant impact, uv_signature (Dipyrimidine transversion point mutations))

=item B<--exclude>

specifies the database filters which the user wants to exclude from his analysis, cannot be used with --include unless --include is set to 'all'.
supported choices are : gene (exclude gene annotation), snp132 (exclude dbSNP_132 filter), snp138 (exclude dbSNP_138 filter), snp138common (exclude common SNPs from dbSNP_138 filter), snp138flagged (exclude flagged SNPs from dbSNP_138 filter), snp142 (exclude dbSNP_142 filter), snp142common (exclude common SNPs from dbSNP_142 filter), snp142flagged (exclude flagged SNPs from dbSNP_142 filter), 1kg (exclude 1000 Genome filter), 1kg_afr (exclude filter of allele frequency of African population of 1000genomes project), 1kg_amr (exclude filter of allele frequency of Ad Mixed American population of 1000genomes project), 1kg_eur (exclude filter of allele frequency of European population of 1000genomes project), 1kg_eas (exclude filter of allele frequency of East Asian population of 1000genomes project), 1kg_sas (exclude filter of allele frequency of South Asian population of 1000genomes project), exac (exclude filter), exac_afr (exclude filter), exac_amr (exclude filter), exac_eas (exclude filter), exac_fin (exclude filter), exac_nfe (exclude filter), exac_oth (exclude filter), exac_sas (exclude filter), cg (exclude Complete Genomics filter), gme (exclude Great Middle East frequencies), dbnsfp (exclude dbnsfp33a database filter : excluding whole-exome SIFT scores, PolyPhen2 HDIV scores, PolyPhen2 HVAR scores, LRT scores, MutationTaster scores, MutationAssessor score, FATHMM scores, PROVEAN scores, VEST scores, CADD scores, DANN scores, M-CAP scores, Eigen scores, GenoCanyon scores, Interpro domains, GTex annotations, FATHMMMKL scores, MetaSVM scores, MetaLR scores, INTEGRATED scores, GERP++ scores, PhyloP scores, Phast scores and SiPhy scores from dbnsfp (dbnsfp version 3.3a - 20170221), sift (exclude SIFT filter), pp2 (exclude PolyPhen.v2 HumanDiv & HumanVar filters), mt (exclude Mutation Taster filter), lrt (exclude LRT filter), phylop (exclude PhyloP filter), gerp_gt2 (whole-genome gerp++ scores greater than 2), gerp_elem (conserved genomic regions by gerp++), ma (exclude filter), fathmm (exclude filter), rsvm (exclude filter), lr (exclude filter), vest (exclude filter), cadd (exclude filter), gerp (exclude filter), siphy (exclude filter), esp6500 (exclude esp6500 filter), esp6500aa (exclude esp6500aa filter), esp6500ea (exclude esp6500ea filter), mirna (exclude mirna region annotation), mirnatarget (exclude mirnatarget region annotation), rnaseq_cuff (exclude Broad Cufflinks RNASeq alignment region annotation), rnaseq_scri (exclude Broad Scripture RNASeq alignment region annotation), clinvar (exclude NCBI\'s ClinVar annotation), cosmic (exclude Sanger's COSMIC (Catalogue of Somatic Mutations in Cancer) somatic disease mutations description region annotation), icgc21 (International Cancer Genome Consortium version 21), nci60 (exclude NCI-60 allele frequency filter), gtx_chipseq (exclude GenomeTrax Predicted ChIP-Seq TFBS description region annotation), gtx_cpg (exclude GenomeTrax CpG Islands description region annotation), gtx_disease (exclude GenomeTrax Disease associations description region annotation), gtx_dnase (exclude GenomeTrax Dnase biding site), gtx_drug (exclude GenomeTrax Drug Targets description region annotation), gtx_gwas (exclude GenomeTrax GWAS Catalogue description region annotation), gtx_hgmd (exclude GenomeTrax HGMD inherited (germ-line) disease mutations description region annotation), gtx_hgmd_disgene (exclude GenomeTrax HGMD disease genes description region annotation), gtx_hgmdimputed (exclude GenomeTrax HGMD IMPUTED mutation), gtx_microsat (exclude GenomeTrax Microsatellite repeats description region annotation), gtx_mirna (exclude GenomeTrax microRNA sequences description region annotation), gtx_omim (exclude GenomeTrax NCBI\'s OMIM description region annotation), gtx_orpha(exclude GenomeTrax rare diseases information), gtx_path (exclude GenomeTrax Pathways membership description region annotation), gtx_pgmd (exclude GenomeTrax Pharmacogenomics Variants description region annotation), gtx_ptms (exclude GenomeTrax Post translational modifications description region annotation), gtx_transfac (exclude GenomeTrax TRANSFAC experimentally verified TFBS description region annotation), gtx_tss (exclude GenomeTrax TSSs (transcription start sites) description region annotation), gtx_vista (exclude vista track for predicted damaging effect of a variation occurring in an enhancer site), enc_hmm (exclude Encode Chromatin States (UCSC track) description region annotation), enc_tfbs (exclude Encode Transcription Factor Binding Sites (UCSC track) description region annotation), simple_repeat (exclude Simple Tandem Repeats (possibly imperfect repeats) located by Tandem Repeats Finder (TRF) (UCSC simpleRepeat track) description region annotation), segment_dups (exclude regions detected as putative genomic duplications (UCSC  genomicSuperDups track) description region annotation), rmsk (exclude regions detected by RepeatMasker as interspersed repeats and low complexity DNA sequencesinterspersed repeats and low complexity DNA sequences), gonl (exclude GoNL project frequencies), wgsdann (exclude dann whole-genome variants score), gwava (exclude gwava region/tss score), wgscadd (exclude cadd/cadd_Phred whole genome variants score), wgsfathmm (exclude fathmm whole-genome coding/noncoding score), hrcr1 (exclude 40M variants from 32K samples HRC_[AF,AC,AN] AND HRC_nonKG_[AF,AC,AN]), dbscsnv11 (exclude dbscsnv splice site effect prediction), exac_nontcga (exclude exac nontcga allele frequency), exac_nonpsych (exclude exac nonpsych allele frequency), kaviar (exclude kaviar allele frequency), wgseigen (exclude eigen wgs coding & noncoding score), rvis_exac (exclude RVIS (Residual Variation Intolerance Score) score based on ExAC sequencing data), revel (exclude REVEL score for predicting pathogenicity of missense variants based on a combination of scores from 13 individual tools), intervar (exclude INTERVAR score for clinical interpretation of missense variants), ebi3222pak (exclude MAF from EBI dataset for 3222 healthy pakistani controls), mpc_exac (ExAC regional missense constraint score), uorf (5'UTR variant impact, uv_signature (Dipyrimidine transversion point mutations))


doIt(-elem => "dbnsfp") || =cut
