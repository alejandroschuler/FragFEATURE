excludes = []

# Water
excludes.extend(["HOH","DOD"])

# Unknown ligand or ion or residue-ligand
excludes.extend(["UNL","UNK","UNX"])

# http://www.sigmaaldrich.com/life-science/proteomics/protein-structural-analysis/xray-crystallography/basic-crystallography-kit.html
# Crystallography Buffers
# ACETATE(2), ADA, BICINE, BIS-TRIS, CACODYLATE(2), CITRATE(3), FORMATE, GLYCINE, IMIDAZOLE, HEPES, MES, MOPS, PHOSPHATE(4), TARTRATE(3), TRIS
excludes.extend(["ACT","ACY","MHA","BCN","BTB","CAC","CAD","FLC","ICT","CIT","FMT","GLY","IMD","EPE","MES","MPO","PO4","2HP","IPS","PI","TAR","TLA","SRT","TRS"])

# Crystallography Precipitating Salts 
# CTAB, SULFATE
excludes.extend(["16A","SO4"])

# Crystallography Precipitating Organic Solvents
# 1,6-HEXANEDIOL, 2-PROPANOL, DIOXANE, DMSO, ETHANOL, ETHYLENE GLYCOL, GLYCEROL, MPD(R/S), TERT-BUTANOL
excludes.extend(["HEZ","IPA","DIO","DMS","EOH","EDO","GOL","MRD","MPD","TBU"])

# Crystallography PEGs & JEFFAMINE
excludes.extend(["1PE","1PG","12P","15P","16P","2PE","6JZ","7PE","7PG","AE3","AE4","ETE","M2M","ME2","P15","P33","P3G","P4G","P6G","PE3","PE4","PE5","PE8"])
excludes.extend(["PEG","PEU","PG0","PG4","PG5","PG6","PGE","PGF","TOE","XPE","JEF"])

# Crystallography Polyamines
# SPERMINE(2), SPERMIDINE 
excludes.extend(["SPK","SPM","SPD"])

# http://www.piercenet.com/browse.cfm?fldID=5558F7E4-5056-8A76-4E55-4F3977738B63
# Common Detergents
# BOG, CHAPS, SDS, TRITON(4)
excludes.extend(["BOG","CPS","SDS","DR6","EGC","OXN","TRT"])

# http://www.science.smith.edu/departments/Biochem/Biochem_353/Common_Buffers.html
# Other Buffers
# BIS-TRIS PROPANE, BORIC ACID, CAPS, CHES, GLYCINE AMIDE, H2CO3/NaHCO3, HEPPSO, PIPES, TAPS, TES, TRIETHANOLAMINE
excludes.extend(["B3P","BO3","CXS","NHE","GM1","CO3","BCT","250","PIN","T3A","NES","211"])

# Reducing Agents
# BETA-MERCAPTO, DTT, TCEP
excludes.extend(["BME","DTT","DTU","DTV","D1D","DTD","TCE"])

# Single heavy atom ligands
excludes.extend(["AG" ,"AL" ,"AR" ,"ARS","AU" ,"AU3","BA" ,"BR" ,"CA" ,"CD" ,"CE" ,"CH2",
                 "CL" ,"CO" ,"3CO","CR" ,"CS" ,"CU" ,"CU1","CU3","ER3","EU" ,"EU3","F"  ,
                 "FE" ,"FE2","GA" ,"GD" ,"GD3","H2S","HG" ,"HO" ,"HO3","IN" ,"IOD","IR" ,
                 "IR3","K"  ,"KR" ,"LA" ,"LI" ,"LU" ,"MG" ,"MN" ,"MN3","MO" ,"4MO","6MO",
                 "NA" ,"NH2","NH3","NH4","NI" ,"3NI","O"  ,"OH" ,"OS" ,"OS4","PB" ,"PD" ,
                 "PR" ,"PT" ,"PT4","RB" ,"RE" ,"RH" ,"RH3","RU" ,"SB" ,"SE" ,"SM" ,"SR" ,
                 "TB" ,"TE" ,"TL" ,"U1" ,"V"  ,"W"  ,"XE" ,"Y1" ,"YB" ,"YB2","YT3","ZN"])
# Two heavy atom ligands
excludes.extend(["06C","0QE","6WO","BMM","C1O","CF0","CMO","CUA","CYN","EHN","FOR","FS1",
                 "HDZ","HOA","MEE","MH2","MMC","MOH","NME","NO","OFE","OXY","PEO","PER","SX"])

# Acetyls to exclude for acetylation sites
acetyls = ["ACE"]

# Chromophores to exclude for chromophore sites
chromos = ["0YG","5ZA","AYG","C12","C99","CCY","CFY","CH6","CH7","CJO","CLV","CR0","CR2",
           "CR5","CR7","CR8","CRF","CRG","CRK","CRO","CRQ","CRU","CRW","CRX","CSH","CSY",
           "CWR","CZO","DYG","GYC","GYS","IEY","IIC","MDO","MFC","NRQ","NYC","NYG","PIA",
           "QLG","RC7","X9Q","XXY","XYG"]

# Sugars to exclude for glycosylation sites
sugars = ["0AT","0BD","10M","147","149","14T","16G","1AR","1GL","1GN","1NA","289","291",
          "293","2DG","2FG","2FL","2FP","2GL","2M4","2M5","3CM","3FM","3LR","3MF","3MG",
          "3SA","46D","46M","4GC","5RP","5SA","5SP","6PG","6SA","7JZ","A2G","AAB","AAL",
          "AAO","ABC","ABD","ABE","ABF","ACG","ACI","ACR","ACX","ADA","ADG","ADR","AF1",
          "AFL","AFO","AFP","AGC","AGH","AGL","AHR","AIG","ALL","AMU","AOG","AOS","ARA",
          "ARB","ARE","ARI","ASG","ASO","AXP","AXR","B0D","B2G","B4G","B6D","B8D","B9D",
          "BCD","BDG","BDP","BDR","BEM","BFP","BGC","BGL","BGP","BGS","BHG","BMA","BNG",
          "BOG","BRI","BXP","BXX","BXY","C3X","C4X","C5X","CAP","CBI","CBS","CDR","CEG",
          "CGF","CHO","CR1","CR6","CRA","CTO","CTR","CTT","D6G","DAF","DAG","DDA","DDB",
          "DDL","DEL","DFR","DFX","DG0","DGC","DGD","DGM","DGS","DIG","DLF","DLG","DMU",
          "DNO","DOM","DP5","DQQ","DQR","DR2","DR3","DR4","DRI","DSR","DT6","E4P","EAG",
          "EBG","EBQ","EJT","EPG","ERE","ERI","F1P","F1X","F6P","FBP","FCA","FCB","FCT",
          "FDP","FDQ","FFC","FIX","FRU","FU4","FUB","FUC","FUD","FUL","FXP","G16","G1P",
          "G2F","G3I","G4D","G4S","G6D","G6P","G6S","GAC","GAD","GAL","GC1","GC4","GCD",
          "GCN","GCO","GCS","GCT","GCU","GCV","GCW","GE1","GFG","GFP","GL0","GL2","GL5",
          "GL6","GL7","GL9","GLA","GLB","GLC","GLD","GLF","GLG","GLO","GLP","GLS","GLT",
          "GLW","GMH","GN1","GP1","GP4","GPH","GPM","GS1","GS4","GSA","GSD","GTE","GTH",
          "GTR","GTZ","GU0","GU1","GU2","GU3","GU4","GU5","GU6","GU8","GU9","GUP","GYP",
          "GYV","H2P","HDL","HMS","HSD","HSG","HSH","HSQ","HSR","HSU","HSX","HSY","HSZ",
          "IAB","IDG","IDR","IDS","IDT","IDU","IDX","IMK","IN1","IPT","KBG","KDA","KDM",
          "KDO","KFN","KO1","KO2","L6S","LAG","LAI","LAK","LAT","LB2","LBT","LDY","LGC",
          "LGU","LM2","LMT","LMU","LOX","LXB","LXZ","M3M","M6P","M8C","MA1","MA2","MA3",
          "MAB","MAG","MAL","MAN","MAT","MAV","MAW","MBG","MDA","MDM","MDP","MFA","MFB",
          "MFU","MGA","MGL","MLB","MMA","MMN","MRP","MTT","MUG","MXY","NAA","NAG","NBG",
          "NDG","NED","NG1","NG6","NGA","NGC","NGE","NGL","NGS","NGY","NHF","NTF","NTO",
          "NTP","NXD","NYT","OPG","OPM","ORP","OX2","P3M","P6P","PA5","PNA","PNG","PNW",
          "PRP","PSJ","PSV","QDK","QPS","QV4","R1P","R1X","R2B","RAA","RAE","RAF","RAM",
          "RAO","RAT","RDP","REL","RER","RG1","RGG","RHA","RIB","RIP","RNS","RNT","ROB",
          "RPA","RUB","S6P","SCR","SGA","SGC","SGN","SGS","SHB","SIO","SOE","SOL","SSG",
          "SUC","SUP","SUS","T6P","TAG","TCB","TDG","TH1","TIA","TM5","TM6","TMR","TMX",
          "TOA","TOC","TRE","TYV","UCD","VG1","X2F","X4S","X5S","XBP","XDN","XDP","XIF",
          "XIM","XLF","XLS","XMM","XUL","XYP","XYS"]
