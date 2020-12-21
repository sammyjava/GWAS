package org.ncgr.gwas;

/**
 * Simple class representing a genomic region.
 */
public class Region {
    String contig;
    int start;
    int end;

    public Region(String region) {
        String[] parts1 = region.split(":");
        this.contig = parts1[0];
        String[] parts2 = parts1[1].split("-");
        this.start = Integer.parseInt(parts2[0]);
        this.end = Integer.parseInt(parts2[1]);
    }

    public Region(String contig, int start, int end) {
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

    public String toString() {
        return contig+":"+start+"-"+end;
    }
}
