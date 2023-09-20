//
// This file holds several Groovy functions
//
@Grab(group='org.codehaus.groovy', module='groovy-yaml', version='3.0.16')
import groovy.yaml.YamlSlurper


class PooledUtils {
    /**
    * Get a List of all the keys in the sample meta yml file 
    * that are marked as 'is_ref'.
    */
    public static List findReferences(sample_meta){
        def samples = new YamlSlurper().parse(sample_meta)["samples"]
        
        List l = samples.findAll {it.value.is_ref}
                .collect{entry -> entry.key}
        return l
        
    }
    public static List createMageckTestArgumentListIllumina(sample_meta){
        def samples = new YamlSlurper().parse(sample_meta)
        
        def mapByGroup= samples.collectEntries{key, value-> [key, [group: value.group, reference: value.reference] ]}
                        .groupBy {it.value.group}

        def retList=[]
        for(entry in mapByGroup){
            
            def referenceCombinations = [:]
            for(sampleEntry in entry.value){
                
                //if a reference value is null then it is skipped or is the control for that sample
                //so ignore those
                referenceCombinations.putAll(sampleEntry.value.reference.findAll{it.value !=null})
            }
            
            for(refEntry in referenceCombinations){
                println(refEntry.key)
                //we are here if we found non-null reference combinations
                //if all combinations were null or empty then we skip this group 
                retList.add([
                    "group" : entry.key,
                    "run_mageck_test" : true,
                    "c" : mapByGroup[refEntry.value].keySet(),
                    "t": entry.value.keySet(),
                    "controlGroup": refEntry.value,
                    "contrast_prefix": refEntry.key
                ])
            }
        }

        return retList
        
        
    }
    
     /**
    * Create a Map of the arguments for the call to MAGECK test.
    The fields are:
    <ul>
    <li>c - the reference/control columns</li>
    <li>t - the test columns.  Note this is a list, one element per invocation</li>
    </ul>
    */
    public static Map createMageckTestArgumentMap(sample_meta){
        def samples = new YamlSlurper().parse(sample_meta)

        def argumentMap = [:]
        argumentMap["group"] = []
        argumentMap["run_mageck_test"] = false
        def mapByGroup= samples.collectEntries{key, value-> [key, [group: value.group, is_ref: value.is_ref] ]}
                        .groupBy {it.value.group}

        for(entry in mapByGroup){
            def key = entry.key
            
            //groups should either all be references or none should be references
            //so checking for the presence or absence should work
            if(entry.value.any(v -> v.value.is_ref)){
                //this is a reference (control group)
                argumentMap["c"] = entry.value.collect(v-> v.key).join(",")
                argumentMap["controlGroup"] = key
                argumentMap["run_mageck_test"] = true
            }else{
                if(!argumentMap["t"]){
                    argumentMap["t"] = [];
                }
                argumentMap["t"].add( [entry.value.collect(v -> v.key).join(",")])
                argumentMap["group"].add(key)
            }        
        }

        return argumentMap
       
    
    }
    // check if the yaml file has multiple read entries per sample
    public static Boolean multipleRead(sample_meta) {
        def samples = new YamlSlurper().parse(sample_meta)["samples"]
        def multiple_read = false
        samples.each {
            def label = it*.key
            def value = it*.value
            if (samples[label[0]].containsKey("read1") && samples[label[0]].containsKey("read2")) {
                multiple_read = true
            }
        }
        return multiple_read
    }

    public static List getSampleLabelFastqMapForMageck(sample_meta) {
        def samples = new YamlSlurper().parse(sample_meta)["samples"]

        List sample_lable_file_list = [[],[],[],[],[]]

        samples.each {

            def label = it*.key
            def value = it*.value
            if (samples[label[0]].containsKey("read1")) {
                sample_lable_file_list[0].add(label[0])
                // The replicates of the same sample are comma separated,
                // the fastq files are space separated
                def filename_string = ""
                for (file in samples[label[0]].read1) {
                    // add comma to the end of file.toString() to make mageck happy
                    filename_string = filename_string + file.toString()+','
                }
                // remove the last comma
                filename_string = filename_string.substring(0, filename_string.length() - 1)
                sample_lable_file_list[1].add(filename_string)
                sample_lable_file_list[2].add(samples[label[0]].group)
                if (samples[label[0]].containsKey("is_ref")) {
                    sample_lable_file_list[3].add(samples[label[0]].is_ref)
                }                 
            }
            if (samples[label[0]].containsKey("read2")) {
                sample_lable_file_list[0].add(label[0])
                def filename_string = ""
                for (file in samples[label[0]].read2) {
                    filename_string = filename_string + file.toString()+','
                }
                // remove the last comma
                filename_string = filename_string.substring(0, filename_string.length() - 1)
                sample_lable_file_list[1].add(filename_string)

            }
        }

        return sample_lable_file_list
    }

    public static List getCountByRepFilesLabel(out_file) {
        def files = new YamlSlurper().parse(out_file)
        List label_file_list = [[],[]]

        files.each {

            def label = it*.key
            def value = it*.value

            label_file_list[0].add(label.join(''))
            label_file_list[1].add(value.join(''))

        }
        return label_file_list
    }

}

