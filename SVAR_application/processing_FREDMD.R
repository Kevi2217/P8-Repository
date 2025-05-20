# Read in data ------------------------------------------------------------
raw_data<-read.csv(file= paste0(getwd(), "/current.csv"),encoding = "UTF-8") #the file "current.csv" was obtained from https://research.stlouisfed.org/econ/mccracken/fred-databases/ on 9/9/2021
colnames(raw_data)[colnames(raw_data)=="S.P.500"]<-"S&P 500" #renamed S.P.500 -> S&P 500, read.csv() changed it
colnames(raw_data)[colnames(raw_data)=="S.P..indust"]<-"S&P: indust" #renamed S.P..indust -> S&P: indust, read.csv() changed it
colnames(raw_data)[colnames(raw_data)=="S.P.div.yield"]<-"S&P div yield" #renamed S.P.div.yield ->S&P div yield, read.csv() changed it
colnames(raw_data)[colnames(raw_data)=="S.P.PE.ratio"]<-"S&P PE ratio" #renamed S.P.PE.ratio -> S&P PE ratio, read.csv() changed it

# Transformation codes of Bernanke ----------------------------------------

# some of these variables are no longer in the FRED-MD database
variables<-c("IPFPNSS","IPFINAL","IPCONGD","IPDCONGD","IPNCONGD",
             "IPBUSEQ","IPMAT", "IPDMAT","IPNMAT","IPMANSICS",
             "IPB51222S","INDPRO","CUMFNS","NAPM","NAPMPI",
             "RPI","W875RX1","HWI","HWIURATIO","CLF16OV",
             "CE16OV","UNRATE","UEMPMEAN","UEMPLT5","UEMP5TO14",
             "UEMP15OV","UEMP15T26","PAYEMS", "CEU0500000001","USGOOD",
             "CES1021000001","USCONS","MANEMP","DMANEMP","NDMANEMP",
             "SRVPRD","USTPU","USWTRADE","USFIRE","USSERV",
             "USGOVT","AWHMAN","AWOTMAN","NAPMEI","DPCERA3M086SBEA",
             "DDURRG3M086SBEA","DNDGRG3M086SBEA","DSERRG3M086SBEA","HOUST","HOUSTNE",
             "HOUSTMW","HOUSTS","HOUSTW","NAPMII","NAPMNOI",
             "NAPMSDI","A0M008","A0M027","S&P 500","S&P: indust",
             "S&P div yield","S&P PE ratio","EXSZUSx","EXJPUSx","EXUSUKx",
             "EXCAUSx","FEDFUNDS","TB3MS","TB6MS","GS1",
             "GS5","GS10","AAA","BAA","TB3SMFFM",
             "TB6SMFFM","T1YFFM","T5YFFM","T10YFFM","AAAFFM",
             "BAAFFM","M1SL","M2SL","M3SL","M2REAL",
             "BOGMBASE","TOTRESNS","NONBORRES","BUSLOANS","NONREVSL",
             "NAPMPRI","WPSFD49207","WPSFD49502","WPSID61","WPSID62",
             "A0M099","CPIAUCSL","CPIAPPSL","CPITRNSL","CPIMEDSL",
             "CUSR0000SAC","CUSR0000SAD","CUSR0000SAS","CPIULFSL","CUSR0000SA0L2",
             "CUSR0000SA0L5","CES2000000008","CES3000000008","UMCSENTx"
)
DRI_mcgraw_names<-c("IPP", "IPF", "IPC", "IPCD", "IPCN",
                    "IPE", "IPM", "IPMD", "IPMND", "IPMFG",
                    "IPUT", "IP", "IPXMCA", "PMI", "PMP",
                    "GMPYQ ", "GMYXPQ","LHEL", "LHELX", "LHEM",
                    "LHNAG", "LHUR", "LHU680", "LHU5", "LHU14",
                    "LHU15", "LHU26", "LPNAG", "LP", "LPGD",
                    "LPMI", "LPCC", "LPEM", "LPED", "LPEN",
                    "LPSP", "LPTU", "LPT", "LPFR", "LPS",
                    "LPGOV", "LPHRM", "LPMOSA", "PMEMP", "GMCQ",
                    "GMCDQ", "GMCNQ", "GMCSQ", "HSFR", "HSNE",
                    "HSMW", "HSSOU", "HSWST", "PMNV", "PMNO",
                    "PMDEL", "MOCMQ", "MSONDQ", "FSPCOM", "FSPIN",
                    "FSDXP", "FSPXE", "EXRSW", "EXRJAN", "EXRUK",
                    "EXRCAN", "FYFF", "FYGM3", "FYGM6", "FYGT1",
                    "FYGT5", "FYGT10", "FYAAAC", "FYBAAC", "SFYGM3",
                    "SFYGM6", "SFYGT1", "SFYGT5", "SFYGT10", "SFYAAAC",
                    "SFYBAAC", "FM1", "FM2", "FM3", "FM2DQ",
                    "FMFBA", "FMRRA", "FMRNBA", "FCLNQ", "CCINRV",
                    "PMCP", "PWFSA", "PWFCSA", "PWIMSA", "PWCMSA",
                    "PSM99Q", "PUNEW", "PU83", "PU84", "PU85",
                    "PUC", "PUCD", "PUS", "PUXF", "PUXHS",
                    "PUXM", "LEHCC", "LEHM", "HHSNTN")
codes<-c(5,5,5,5,5,
         5,5,5,5,5,
         5,5,1,1,1,
         5,5,5,4,5,
         5,1,1,1,1,
         1,1,5,5,5,
         5,5,5,5,5,
         5,5,5,5,5,
         5,1,1,1,5,
         5,5,5,4,4,
         4,4,4,1,1,
         1,5,5,5,5,
         1,1,5,5,5,
         5,1,1,1,1,
         1,1,1,1,1,
         1,1,1,1,1,
         1,5,5,5,5,
         5,5,5,5,5,
         1,5,5,5,5,
         5,5,5,5,5,
         5,5,5,5,5,
         5,5,5,1
)
bernanke_codes<-data.frame(cbind(variables,codes, DRI_mcgraw_names))

# Define the slow and fast variables --------------------------------------

# slow & fast variables matched to the original Bernanke paper
slow_variables<-c("IPFPNSS", "IPFINAL", "IPCONGD", "IPDCONGD", "IPNCONGD", "IPBUSEQ",
                  "IPMAT", "IPDMAT", "IPNMAT", "IPMANSICS", "IPB51222S",
                  #"IPMINE" does't exist
                  "INDPRO", "CUMFNS",
                  #"NAPM", removed from FRED-MD
                  #"NAPMPI", removed from FRED-MD
                  "RPI", "W875RX1", "HWI",
                  "HWIURATIO", "CLF16OV", "CE16OV", "UNRATE", "UEMPMEAN",
                  "UEMPLT5", "UEMP5TO14", "UEMP15OV", "UEMP15T26", "PAYEMS",
                  #"CEU0500000001", doesn't exist
                  "USGOOD", "CES1021000001", "USCONS", "MANEMP", "DMANEMP",
                  "NDMANEMP", "SRVPRD", "USTPU", "USWTRADE", "USFIRE",
                  #"USSERV", doesn't exist
                  "USGOVT", "AWHMAN", "AWOTMAN",
                  #"NAPMEI", removed from FRED-MD
                  "DPCERA3M086SBEA",
                  "DDURRG3M086SBEA", #"DDURRA3M086SBEA", wrong name
                  "DNDGRG3M086SBEA",  #"DNDGRA3M086SBEA", wrong name
                  "DSERRG3M086SBEA", #"DSERRA3M086SBEA", wrong name
                  "WPSFD49207", "WPSFD49502",
                  "WPSID61", "WPSID62",
                  #"A0M099", doesn't exist
                  "CPIAUCSL", "CPIAPPSL", "CPITRNSL",
                  "CPIMEDSL", "CUSR0000SAC", "CUSR0000SAD", "CUSR0000SAS", "CPIULFSL",
                  "CUSR0000SA0L2", "CUSR0000SA0L5", "CES2000000008", "CES3000000008"
)
fast_variables<-c("HOUST", "HOUSTNE", "HOUSTMW", "HOUSTS", "HOUSTW",
                  #"NAPMII", removed from FRED-MD
                  #"NAPMNOI", removed from FRED-Md
                  #"NAPMSDI", removed from FRED-MD
                  #"A0M008", doesn't exist
                  #"A0M027", doesn't exist
                  "S&P 500", "S&P: indust", "S&P div yield",
                  "S&P PE ratio", "EXSZUSx", "EXJPUSx", "EXUSUKx", "EXCAUSx",
                  #"FEDFUNDS" this is the shock variable, kept separately
                  "TB3MS", "TB6MS", "GS1", "GS5", "GS10", "AAA", "BAA", "TB3SMFFM",
                  "TB6SMFFM",
                  "T1YFFM", "T5YFFM", "T10YFFM", "AAAFFM", "BAAFFM",
                  "M1SL", "M2SL",
                  #"M3SL", discontinued
                  "M2REAL",
                  "BOGMBASE", #previously called "AMBSL"
                  "TOTRESNS", "NONBORRES",
                  "BUSLOANS", "NONREVSL"
                  #"NAPMPRI" removed from FRED-MD
                  #"UMCSENTx" #too many missing(230)
)
# new variables from FRED-MD which we decided should be slow, based on which type of variable they are
# General rules following Bernanke:
# Prices: slow (except the NAPM COMMODITY PRICES INDEX by bernanke, and we also chose to exclude OILPRICEx, and have that in fast instead)
# Output & Income: slow
# Labor Market: slow
# Consumption: slow

# Interest & Exchange Rates: fast
# Money & Credit: fast
# Stock Market: fast
# Housing: fast
new_slow_variables<-c("IPFUELS", "UEMP27OV", "CLAIMSx", "USTRADE",
                      "CES0600000007",
                      "PPICMM",
                      "PCEPI",
                      "CES0600000008"
)
new_fast_variables<-c("CMRMTSPLx", "RETAILx",
                      "PERMIT", "PERMITNE", "PERMITMW",
                      "PERMITS",
                      "PERMITW",
                      #"ACOGNO", #too many missing(400)
                      "AMDMNOx",
                      #"ANDENOx",#too many missing(111)
                      "AMDMUOx", "BUSINVx" , "ISRATIOx",
                      "REALLN", "CONSPI",
                      #"MZMSL", discontinued
                      "CP3Mx", "COMPAPFFx",
                      #"TWEXAFEGSMTHx", #previously called "TWEXMMTH", too many missing (170)
                      "DTCOLNVHFNM", "DTCTHFNM", "INVEST",
                      "OILPRICEx"
)
federal_funds<-"FEDFUNDS"

removed_variables<-c("VXOCLSx", #too many missing (43)
                     "UMCSENTx", #too many missing(230)
                     "ACOGNO", #too many missing(400)
                     "ANDENOx", #too many missing(111)
                     "TWEXAFEGSMTHx" # too many missing (170)
)
# Clean the data ----------------------------------------------------------

# filling in the bernanke transformation codes wherever I can match their data to FRED-MD
fred_tcode<-as.numeric(c(1,raw_data[1,-1]))
fred_names<-colnames(raw_data)
b_tcode<-fred_tcode
for(i in 2:length(fred_names)){
  b_id<-which(group_which(fred_names[i],bernanke_codes$variables)==T)
  if(length(b_id)==1){
    if(!(b_tcode[i]==7 && as.numeric(bernanke_codes$codes[b_id])==5)){
      b_tcode[i]<-as.numeric(bernanke_codes$codes[b_id])
    }
  }
}
# to see the differences between the "regular" FRED-MD transformations and what we do 
data.frame(names=fred_names,fred=fred_tcode,bernanke=b_tcode)

# this should give a warning message about deleting 12 rows of missing data - intended behavior
CDbmedium<-clean_data(raw_data=raw_data,
                      slow_names = c(slow_variables, new_slow_variables),
                      FFR_name = "FEDFUNDS",
                      fast_names = c(fast_variables, new_fast_variables),
                      start_date = "1/1/1959",
                      end_date = "10/1/2008",
                      transform_codes = b_tcode
)
saveRDS(CDbmedium, "CDbmedium.RData")
