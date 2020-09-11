package me.songbx.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.readers.TabixReader;

public class Tests2 {

	
	
	public static void main2(String[] args) {
    	try {
    		String bwFilePosition =   "/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/mapreadsToCns/B73.bw";
    		Path bwFile = Paths.get(bwFilePosition);
    		WigFileReader wig = WigFileReader.autodetect(bwFile);
    		for( int position=1; position<88; ++position   ) {
    			Contig result = wig.query("3:29670514-29670601", position, position);
    			double thisMean = result.mean();
    			System.out.println(position + " " + thisMean);
    		}
    		
    		wig.close();
    		
		} catch (IOException e) {
		} catch (WigFileException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		String mer = "ACATAATCAGTTTTACCAGT";
		String rmer = "ACTGGTAAAACTGATTATGT";
		if( mer.compareTo(rmer) < 0 ){
			System.out.println("yes");
		}else {
			System.out.println("no");
		}
	}
}
