//
// This file holds several Groovy functions
//
@Grab(group='org.codehaus.groovy', module='groovy-yaml', version='3.0.15')
import groovy.yaml.YamlSlurper


class PooledUtils {

    // check if the yaml file has multiple read entries per sample
    public static Boolean multipleRead(sample_meta) {
        def samples = new YamlSlurper().parse(sample_meta)
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
        def samples = new YamlSlurper().parse(sample_meta)
        List sample_lable_file_list = [[],[]]

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

