/*
 * Install Interproscan
 */
process installInterproscan {

    input:


    output:
    path("${tarball}")

    script:
    interproscan_link=params.interproscan_link
    tarball=params.interproscan_tarball
    """
    wget ${interproscan_link}
    wget ${interproscan_link}.md5
    md5sum -c ${tarball}.md5
    """
}
