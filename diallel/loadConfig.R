if (configfile.has(config, "datafile")) {datafile <- configfile.string(config, "datafile")}
if (configfile.has(config, "MIMPlist")) {MIMPlist <- configfile.string(config, "MIMPlist")}
if (configfile.has(config, "imputefile")) {imputefile <- configfile.string(config, "imputefile")}
if (configfile.has(config, "reps")) {reps <- configfile.integer(config, "reps")}
if (configfile.has(config, "savedir")) {savedir <- configfile.get(config, "savedir")}
if (configfile.has(config, "strategy")) {strategy <- configfile.strings(config, "strategy")}
if (configfile.has(config, "trt_colname")) {trt_colname <- configfile.get(config, "trt_colname")}
if (configfile.has(config, "trt_string")) {trt_string <- configfile.get(config, "trt_string")}
if (configfile.has(config, "ctrl_string")) {ctrl_string <- configfile.get(config, "ctrl_string")}
if (configfile.has(config, "subset")) {subset <- configfile.get(config, "subset")}
if (configfile.has(config, "subset_colname")) {subset_colname <- configfile.get(config, "subset_colname")}
if (configfile.has(config, "assignment")) {assignment <- configfile.get(config, "assignment")}
if (configfile.has(config, "phenotype")) {phenotype <- configfile.get(config, "phenotype")}
if (configfile.has(config, "phenImpute")) {phenImpute <- configfile.get(config, "phenImpute")}
if (configfile.has(config, "phenbatch")) {phenbatch <- configfile.strings(config, "phenbatch")}
if (configfile.has(config, "phenrandom")) {phenrandom <- configfile.strings(config, "phenrandom")}
if (configfile.has(config, "phenfixed")) {phenfixed <- configfile.strings(config, "phenfixed")}
if (configfile.has(config, "phenType")) {phenType <- configfile.strings(config, "phenType")}
if (configfile.has(config, "phenBS")) {phenBS <- configfile.strings(config, "phenBS")}
if (configfile.has(config, "covImpute")) {covImpute <- configfile.get(config, "covImpute")}
if (configfile.has(config, "covbatch")) {covbatch <- configfile.strings(config, "covbatch")}
if (configfile.has(config, "covrandom")) {covrandom <- configfile.strings(config, "covrandom")}
if (configfile.has(config, "covfixed")) {covfixed <- configfile.strings(config, "covfixed")}
if (configfile.has(config, "covType")) {covType <- configfile.strings(config, "covType")}
if (configfile.has(config, "covBS")) {covBS <- configfile.strings(config, "covBS")}
if (configfile.has(config, "numChains")) {numChains <- configfile.integer(config, "numChains")}
if (configfile.has(config, "lengthChains")) {lengthChains <- configfile.integer(config, "lengthChains")}
if (configfile.has(config, "burnin")) {burnin <- configfile.integer(config, "burnin")}
if (configfile.has(config, "thin")) {thin <- configfile.integer(config, "thin")}
if (configfile.has(config, "contrast")) {contrast <- configfile.get(config, "contrast")}
if (configfile.has(config, "type")) {type <- configfile.get(config, "type")}
if (configfile.has(config, "batch")) {batch <- configfile.strings(config, "batch")}
if (configfile.has(config, "random")) {random <- configfile.strings(config, "random")}
if (configfile.has(config, "fixed")) {fixed <- configfile.strings(config, "fixed")}
if (configfile.has(config, "BS")) {BS <- configfile.get(config, "BS")}
if (configfile.has(config, "num_piles")) {num_piles <- configfile.integer(config, "num_piles")}
