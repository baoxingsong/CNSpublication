package me.songbx.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.readers.TabixReader;

public class Tests {

	public static void main1(String[] args) {
		SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
		String chr = "1";
		int start = 204418759;
		int end = 204518758;
		
		SAMRecordIterator it;
		
		it = reader.queryOverlapping(chr, start, end);
		while (it.hasNext()) {
            final SAMRecord samRecord = it.next();
            System.out.println(samRecord.getReadName() +  " " + samRecord.getCigarString() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd() + " " + samRecord.getReadNegativeStrandFlag());
        }
		
		
		reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
		
		chr = "4";
		start = 51137609;
		end = 51237608;
		it = reader.queryOverlapping(chr, start, end);
		while (it.hasNext()) {
            final SAMRecord samRecord = it.next();
            System.out.println(samRecord.getReadName() +  " " + samRecord.getCigarString() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd() + " " + samRecord.getReadNegativeStrandFlag());
        }
	}
	
	public static void main2(String[] args) {
    	try {
			TabixReader tr = new TabixReader("/media/bs674/panAndAssemblyfi/maize282Genotyping/callablePipeline/Tzi11_callable_status.bed.gz");
			String s;
			String chr = "1";
			int start = 10;
			int end = 10000;
			TabixReader.Iterator iter = tr.query(chr+ ":" + start + "-" + end); // get the iterator
			while ((s = iter.next()) != null) {
				System.out.println(s); //haha
				// here should be tested in two ways, 1) mean value are different; 2) the non-mappable accessions have higher variance than mappable accessions
			}
			tr.close();
		} catch (IOException e) {
		}
	}
	
	
	public static void main(String[] args) {
		int length = 10;
		String s= new String(new char[length]).replace("\0", "-");
		System.out.println(s);
		
	}
	
	
}
