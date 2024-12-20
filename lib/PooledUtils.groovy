//
// This file holds several Groovy functions
//
@Grab(group='org.codehaus.groovy', module='groovy-yaml', version='3.0.16')

import java.nio.file.Path

import groovy.yaml.YamlSlurper
import nextflow.Nextflow


class PooledUtils {
    /**
    * Validate that a file path exists and return as a Path object,
    * or return an appropriate Nextflow error message.
    */
    public static Path validatePathParam(path, param_name, is_dir = false) {
        try {
            Path file = Nextflow.file(path, checkIfExists: true)
            if (is_dir && !file.isDirectory()) {
                throw new java.nio.file.NotDirectoryException(path)
            }
            return file
        } catch (java.nio.file.NoSuchFileException e) {
            // return the error message if the path does not exist
            Nextflow.error "The path " + e.getMessage() + " does not exist"
        } catch (java.nio.file.NotDirectoryException e) {
            // return the error message if the path is not a directory
            Nextflow.error "The path " + e.getMessage() + " is not a directory!"
        } catch (java.nio.file.AccessDeniedException e) {
            // return the error message if the path cannot be accessed
            Nextflow.error "The path " + e.getMessage() + " cannot be accessed!"
        } catch (e) {
            // return everything else
            Nextflow.error "The ${param_name} parameter definition failed with following reason " + e
        }
    }

    /**
    * Create a list of Maps containing the values needed to run MAGeCK test for
    * each contrast.
    */
    public static List createMageckTestContrasts(sample_meta){
        def samples = new YamlSlurper().parse(sample_meta)
        
        def samplesByGroup = samples.collectEntries { key, value ->
            [key, [group: value.group, reference: value.reference, is_ref: value.is_ref]]
        }
        .groupBy {it.value.group}

        // In case the samples do not contain a reference key defining the reference
        // groups, find the reference group with the is_ref key (ie, the group that
        // has a sample with is_ref set to 1)
        // TODO: Deprecate the is_ref approach to the reference group in the metadata
        def refGroup = samplesByGroup.findAll{it.value.any{it.value.is_ref}}*.key[0]

        samplesByGroup.collectMany { group, groupSamples ->
            // For this group, get a map of the analyses and corresponding reference groups
            def analysis_refs = groupSamples.collect { _, value ->
                if (!value.reference) {
                    // TODO: Remove this after is_ref usage is deprecated
                    value.reference = ['': (group != refGroup ? refGroup : null)]
                }
                // Analyses with no reference value should be ignored (either this
                // group is not part of that analysis or it is the reference group)
                value.reference.findAll {it.value != null}
            }
            .inject { a, b -> a + b }

            // For each analysis, create a map with the inputs needed for running MAGeCK
            // test to contrast this group to the appropriate reference control group
            analysis_refs.collect { analysis_name, reference ->
                [
                    group: group,
                    samples: groupSamples.keySet(),
                    refGroup: reference,
                    refSamples: samplesByGroup[reference].keySet(),
                    contrast: analysis_name
                ]
            }
        }
    }

    public static Map getSamplesFastqMap(meta_yaml, fastq_dir) {
        def samples = new YamlSlurper().parse(meta_yaml)["samples"]

        samples.collect { sample_name, metadata ->
            def fastqs = [metadata.read1]
            if (metadata.read2) {
                fastqs += [metadata.read2]
            }

            fastqs
            .transpose()
            .collect {
                [
                    meta: [id: sample_name],
                    fastqs_r1: [fastq_dir.resolve(it[0])],
                    fastqs_r2: it[1] ? [fastq_dir.resolve(it[1])] : [],
                    single_end: it.size() == 1
                ]
            }
        }
        .flatten()
        .inject { a, b ->
            // Throw an error if the samples are a mix of single-end and paired-end
            if (a.single_end != b.single_end) {
                error "Samples contain a mix of single-end and paired-end data."
            }
            [
                meta: [id: "${a.meta.id},${b.meta.id}"],
                fastqs_r1: a.fastqs_r1 + b.fastqs_r1,
                fastqs_r2: a.fastqs_r2 + b.fastqs_r2,
                single_end: b.single_end
            ]
        }
    }

    public static List getRepresentationBasenames(repr_basename_yaml) {
        def basenames = new YamlSlurper().parse(repr_basename_yaml)
        basenames.collect { representation, basename -> [basename, representation] }
    }

}
