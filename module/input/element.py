SFED_radi = {
    "1": 1.181,
    "2": 1.052,
    "3": 1.038,
    "4": 1.096,
    "5": 1.08,
    "6": 1.189,
    "7": 1.24,
    "8": 0.9,
    "9": 0.95,
    "11": 1.947,
    "12": 1.886,
    "13": 1.833,
    "14": 1.855,
    "15": 1.817,
    "16": 1.734,
    "17": 1.742,
    "18": 1.801,
    "21": 1.545,
    "22": 1.518,
    "23": 1.518,
    "24": 1.518,
    "25": 1.517,
    "26": 1.429,
    "27": 1.429,
    "31": 1.67,
    "32": 1.66,
    "33": 1.66,
    "34": 1.704,
    "35": 1.696,
    "36": 1.68,
    "37": 1.645,
    "38": 1.495,
    "39": 1.69,
    "41": 3.749,
    "42": 2.101,
    "51": 1.47,
    "52": 1.7,
    "53": 1.89,
    "54": 1.9,
}

# refrelective index, hydrogen bond acidity (donor), basicity (acceptor), surface tension, dielectric constant
solv_prop_nomalized = {
    "1,2,4-trimethylbenzene": [0.922058824, 0, 0.157024793, 0.369667046, 0.01302478],
    "1,2-dibromoethane": [
        0.942830882,
        0.12195122,
        0.140495868,
        0.500696291,
        0.027154736,
    ],
    "1,2-dichloroethane": [
        0.885294118,
        0.12195122,
        0.090909091,
        0.403342195,
        0.055754405,
    ],
    "1,4-dioxane": [0.871568627, 0, 0.52892562, 0.414609444, 0.012169053],
    "1-bromooctane": [0.88995098, 0, 0.099173554, 0.363083935, 0.027667401],
    "1-butanol": [0.857414216, 0.451219512, 0.396694215, 0.315736169, 0.095442181],
    "1-chlorobutane": [0.859068627, 0, 0.082644628, 0.293454868, 0.039046684],
    "1-decanol": [0.880637255, 0.451219512, 0.396694215, 0.360931764, 0.041467511],
    "1-fluorooctane": [0.853860294, 0, 0.082644628, 0.297912397, 0.021420705],
    "1-heptanol": [0.872671569, 0.451219512, 0.396694215, 0.360805165, 0.062340308],
    "1-hexanol": [0.868995098, 0.451219512, 0.396694215, 0.326750222, 0.068888767],
    "1-iodohexadecane": [0.907230392, 0, 0.123966942, 0.408785922, 0.019459251],
    "1-nonanol": [0.878553922, 0.451219512, 0.396694215, 0.353082669, 0.047351872],
    "1-octanol": [0.875612745, 0.451219512, 0.396694215, 0.343081403, 0.054311123],
    "1-pentanol": [0.864031863, 0.451219512, 0.396694215, 0.321053298, 0.083314978],
    "1-propanol": [0.848039216, 0.451219512, 0.396694215, 0.295227244, 0.113015969],
    "2,2,4-trimethylpentane": [0.852634804, 0, 0, 0.213666287, 0.010659692],
    "2,6-dimethylpyridine": [0.916237745, 0, 0.520661157, 0.392065071, 0.039501652],
    "2-butanone": [0.844852941, 0, 0.421487603, 0.303456134, 0.100471916],
    "2-methoxyethanol": [
        0.859313725,
        0.365853659,
        0.694214876,
        0.390429168,
        0.094713656,
    ],
    "2-methylpyridine": [0.916482843, 0, 0.479338843, 0.417774402, 0.054808921],
    "2-pentanol": [0.859681373, 0.402439024, 0.462809917, 0.296873022, 0.072522026],
    "3-methyl-1-butanol": [
        0.862132353,
        0.451219512,
        0.396694215,
        0.300164578,
        0.080947137,
    ],
    "4-methyl-2-pentanone": [0.855514706, 0, 0.421487603, 0.292916825, 0.070964207],
    "acetic_acid": [0.840686275, 0.743902439, 0.363636364, 0.343081403, 0.034431718],
    "acetone": [0.832598039, 0.048780488, 0.404958678, 0.29699962, 0.112848568],
    "acetonitrile": [0.823651961, 0.085365854, 0.26446281, 0.362830738, 0.196520374],
    "aniline": [0.971997549, 0.317073171, 0.338842975, 0.533232055, 0.037930617],
    "anisole": [0.929779412, 0, 0.239669421, 0.444360046, 0.023263767],
    "benzene": [0.919791667, 0, 0.115702479, 0.357260413, 0.012503304],
    "benzonitrile": [0.93682598, 0, 0.272727273, 0.49107482, 0.14092511],
    "benzyl_alcohol": [0.943382353, 0.402439024, 0.462809917, 0.440530447, 0.068595264],
    "bromobenzene": [0.955698529, 0, 0.074380165, 0.446132422, 0.029710352],
    "bromoethane": [0.872487745, 0, 0.099173554, 0.299025193, 0.049614537],
    "bromoform": [0.980698529, 0.182926829, 0.049586777, 0.568046588, 0.023396476],
    "carbon_disulfide": [0.999938725, 0, 0.05785124, 0.399797443, 0.014375],
    "carbon_tetrachloride": [0.894669118, 0, 0, 0.334599316, 0.012268722],
    "chlorobenzene": [0.933884804, 0, 0.05785124, 0.417647804, 0.031370044],
    "chloroform": [0.885968137, 0.182926829, 0.016528926, 0.337637676, 0.025943282],
    "chlorohexane": [0.870036765, 0, 0.082644628, 0.325737435, 0.032759361],
    "cyclohexane": [0.874142157, 0, 0, 0.312064818, 0.011104075],
    "cyclohexanone": [0.888909314, 0, 0.462809917, 0.437650335, 0.086005507],
    "decalin": [0.890196078, 0, 0, 0.39998734, 0.012092511],
    "dichloromethane": [0.872671569, 0.12195122, 0.041322314, 0.344347386, 0.049174009],
    "diethyl_ether": [0.82879902, 0, 0.338842975, 0.210786175, 0.023348018],
    "diisopropyl_ether": [0.83817402, 0, 0.338842975, 0.21863527, 0.018612335],
    "dimethylacetamide": [0.881127451, 0, 0.644628099, 0.419356881, 0.208043502],
    "dimethylformamide": [0.876531863, 0, 0.611570248, 0.457906064, 0.204950441],
    "dimethylsulfoxide": [0.905330882, 0, 0.727272727, 0.543359919, 0.257850771],
    "di-n-butyl_ether": [0.857352941, 0, 0.371900826, 0.284149892, 0.016780286],
    "diphenyl_ether": [0.857843137, 0, 0.165289256, 0.308266869, 0.020539648],
    "ethanol": [0.834007353, 0.451219512, 0.396694215, 0.278136473, 0.13685022],
    "ethoxybenzene": [0.92377451, 0, 0.26446281, 0.410305102, 0.023015969],
    "ethyl_acetate": [0.840870098, 0, 0.371900826, 0.296113432, 0.03296641],
    "ethylbenzene": [0.916605392, 0, 0.123966942, 0.363970123, 0.013402533],
    "ethylene_glycol": [
        0.877328431,
        0.707317073,
        0.644628099,
        0.607545259,
        0.221616189,
    ],
    "fluorobenzene": [0.899754902, 0, 0.082644628, 0.337511077, 0.029845815],
    "iodobenzene": [0.992647059, 0, 0.099173554, 0.490062033, 0.025038546],
    "isobutanol": [0.855085784, 0.451219512, 0.396694215, 0.285384226, 0.092382159],
    "isopropanol": [0.844117647, 0.402439024, 0.462809917, 0.264970249, 0.106082048],
    "isopropylbenzene": [0.913909314, 0, 0.132231405, 0.350487403, 0.013057269],
    "m-cresol": [0.945955882, 0.695121951, 0.280991736, 0.451829345, 0.068502203],
    "mesitylene": [0.91875, 0, 0.157024793, 0.348778326, 0.012472467],
    "methanol": [0.814215686, 0.524390244, 0.388429752, 0.279402456, 0.179587933],
    "methyl_acetate": [0.834191176, 0, 0.371900826, 0.313077605, 0.03778359],
    "methyl_phenyl_ketone": [0.941911765, 0, 0.396694215, 0.494239777, 0.096035242],
    "methyl_tert-butyl_ether": [0.838848039, 0, 0.454545455, 0.224711989, 0.024779736],
    "methylformamide": [
        0.877389706,
        0.487804878,
        0.454545455,
        0.489935435,
        0.999790198,
    ],
    "m-xylene": [0.917401961, 0, 0.132231405, 0.376756551, 0.012928414],
    "n-butyl_acetate": [0.854227941, 0, 0.371900826, 0.314976579, 0.027500551],
    "n-butylbenzene": [0.912867647, 0, 0.123966942, 0.363558678, 0.012995595],
    "n-decane": [0.864093137, 0, 0, 0.284149892, 0.010928414],
    "n-dodecane": [0.871078431, 0, 0, 0.315356374, 0.011046256],
    "n-heptane": [0.850367647, 0, 0, 0.248765667, 0.01052478],
    "n-hexadecane": [0.878982843, 0, 0, 0.342448411, 0.011234581],
    "n-hexane": [0.842463235, 0, 0, 0.226484365, 0.010362885],
    "nitrobenzene": [0.953553922, 0, 0.231404959, 0.578047854, 0.191680066],
    "nitroethane": [0.852757353, 0.024390244, 0.272727273, 0.406760349, 0.155779736],
    "nitromethane": [0.846629902, 0.073170732, 0.256198347, 0.462463603, 0.201334251],
    "n-methyl-2-pyrrolidone": [0.900735294, 0, 0.652892562, 0.515255096, 0.177312775],
    "n-nonane": [0.861151961, 0, 0, 0.283327003, 0.010795705],
    "n-octane": [0.85625, 0, 0, 0.267628814, 0.010686123],
    "n-pentadecane": [0.877144608, 0, 0, 0.337226231, 0.011196586],
    "n-pentane": [0.831801471, 0, 0, 0.195562729, 0.010116189],
    "n-undecane": [0.882230392, 0, 0, 0.306494493, 0.010963656],
    "o-dichlorobenzene": [0.95067402, 0, 0.033057851, 0.453221927, 0.055037996],
    "o-nitrotoluene": [0.946691176, 0, 0.223140496, 0.521205216, 0.14135022],
    "perfluorobenzene": [0.844178922, 0, 0, 0.278516268, 0.011172907],
    "p-isopropyltoluene": [0.913541667, 0, 0.157024793, 0.336755184, 0.01229185],
    "propylene_carbonate": [0.86942402, 0, 0.495867769, 0.530826687, 0.354625551],
    "pyridine": [0.924938725, 0, 0.429752066, 0.462843398, 0.071462555],
    "sec-butanol": [0.856495098, 0.402439024, 0.462809917, 0.297126219, 0.087795154],
    "sec-butylbenzene": [0.912683824, 0, 0.132231405, 0.354886695, 0.012910793],
    "tert-butanol": [0.850367647, 0.37804878, 0.495867769, 0.252690214, 0.068667401],
    "tert-butylbenzene": [0.914644608, 0, 0.132231405, 0.349886062, 0.012911344],
    "tetrachloroethene": [0.922365196, 0, 0, 0.402582605, 0.012488987],
    "tetrahydrofuran": [0.860906863, 0, 0.396694215, 0.335485504, 0.040890419],
    "tetrahydrothiophene-1,1-dioxide": [
        0.908884804,
        0,
        0.727272727,
        0.449423978,
        0.242082599,
    ],
    "tetralin": [0.94442402, 0, 0.157024793, 0.419340133, 0.015258811],
    "toluene": [0.916727941, 0, 0.115702479, 0.353589062, 0.013073238],
    "tributylphosphate": [0.871568627, 0, 1, 0.241982403, 0.04503359],
    "triethylamine": [0.858455882, 0, 0.652892562, 0.25598177, 0.013123348],
    "water": [0.81495098, 1, 0.289256198, 0.911381187, 0.431471916],
    "logp": [
        -0.044473435,
        0.402331736,
        -0.078766874,
        0.416642071,
        0.276510845,
    ],  # (water-1-octanol)/1.364
    "pampa": [
        -0.044473435,
        0.402331736,
        -0.078766874,
        0.416642071,
        0.276510845,
    ],  # (water-1-octanol)/1.364
}
