/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER HELPERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    On Nextflow >= 26.04 a CLI parameter arrives as a String, even an explicit
    boolean: `--include_binning false` is the String "false", which is truthy
    in Groovy, so `if (params.include_binning)` would run the branch
    (Nextflow 25.10 still converts it to a Boolean). nf-schema validates the
    value but does not convert it, and a CLI value overrides any conversion
    written in nextflow.config, so every boolean parameter is read through
    flagOn() at the point of use (design doc Q13). Modules and config files,
    which cannot include this file, use the same expression inline:
    params.<name>.toString().toBoolean()
----------------------------------------------------------------------------------------
*/

// true for Boolean true and the String "true"; false for false, "false" and
// null (nf-schema has already rejected any other value)
def flagOn(value) {
    return value.toString().toBoolean()
}
