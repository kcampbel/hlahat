/*
 * -----------------------------------------------------
 *  Utility functions 
 * -----------------------------------------------------
 */

/*
 * Extract name of software tool from process name using $task.process
 */
def getSoftwareName(task_process) {
    return task_process.tokenize(':')[-1]
}

/*
 * Function to initialise default values and to generate a Groovy Map of available options for nf-core modules
 */
def initOptions(Map args) {
    def Map options = [:]
    options.args          = args.args ?: ''
    options.args2         = args.args2 ?: ''
    options.publish_by_id = args.publish_by_id ?: false
    options.publish_dir   = args.publish_dir ?: ''
    options.publish_files = args.publish_files
    options.suffix        = args.suffix ?: ''
    return options
}

/*
 * Tidy up and join elements of a list to return a path string
 */
def getPathFromList(path_list) {
    def paths = path_list.findAll { item -> !item?.trim().isEmpty() }  // Remove empty entries
    paths = paths.collect { it.trim().replaceAll("^[/]+|[/]+\$", "") } // Trim whitespace and trailing slashes
    return paths.join('/')
}

/*
 * Function to save/publish module results
 */
def saveFiles(Map args) {
    if (!args.filename.endsWith('.version.txt')) {
        def ioptions = initOptions(args.options)
        def path_list = [ ioptions.publish_dir ?: args.publish_dir ]
        if (ioptions.publish_by_id) {
            path_list.add(args.publish_id)
        }
        if (ioptions.publish_files instanceof Map) {
            for (ext in ioptions.publish_files) {
                if (args.filename.endsWith(ext.key)) {
                    def ext_list = path_list.collect()
                    ext_list.add(ext.value)
                    return "${getPathFromList(ext_list)}/$args.filename"
                }
            }
        } else if (ioptions.publish_files == null) {
            return "${getPathFromList(path_list)}/$args.filename"
        }
    }
}

def create_hlacn_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id              = row.sample
    meta.specimen_id     = row.pact_id
    meta.patient_id      = row.patient_info_patient_id
    meta.study           = row.patient_info_study_id
    meta.dob             = row.patient_info_dob
    meta.tcga_study_code = row.patient_info_patient_tumorType

    def data = [:]
    data.normal_bam = file(row.normal_bam)
    data.normal_bai = file(row.normal_bai)
    data.tumor_bam  = file(row.tumor_bam)
    data.tumor_bai  = file(row.tumor_bai)
    data.sequenzaModelRData = file(row.sequenzaModelRData)
    data.epic_hlatypes = file(row.epic_hlatypes)

    array = [ meta, data ]
    return array
}
