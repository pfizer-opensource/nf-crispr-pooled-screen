class OutputFiles {
    // Generate a list of expected output files from running MAGeCK test
    static List<String> mageckTest(String prefix,
                                   Map<String, String> analyses,
                                   List<Integer> representations,
                                   List<String> suffixes) {
        analyses
            .collectMany { analysis, contrasts ->
                [suffixes, representations, contrasts].combinations { suffix, repr, contrast ->
                    [analysis: analysis, contrast: contrast, suffix: suffix, repr: repr]
                }
            }
            .collect {
                def dotAnalysis = it.analysis ? ".${it.analysis}" : ""
                "analysis_${it.repr}X/mageck_test${dotAnalysis}/${prefix}.${it.contrast}.${it.repr}X.${it.suffix}"
            }
    }
}