package me.songbx.impl;

import java.io.File;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.model.Alignment;

public class CigarToAlignment {
	public static Alignment getAlignment( HashMap<String, String> referenceGenome, SAMRecord samRecord ) {
		String referenceSeq = referenceGenome.get(samRecord.getContig()).substring(samRecord.getAlignmentStart()-1, samRecord.getAlignmentEnd());
		String querySeq = samRecord.getReadString();
		String refAln = "";
		String queAln = "";
		int refPosition = 0;
		int queryPostion = 0;
		Pattern pattern = Pattern.compile("(\\d+)(\\w)");
		
		Matcher matcher = pattern.matcher(samRecord.getCigarString());
    	while( matcher.find() ) {
    		if( matcher.group(2).compareTo("H") == 0 || matcher.group(2).compareTo("S") == 0 ) {
    			
    		}else if( matcher.group(2).compareTo("M") == 0 || matcher.group(2).compareTo("=") == 0 || matcher.group(2).compareTo("X") == 0 ) {
    			int length = Integer.parseInt(matcher.group(1));
    			refAln += referenceSeq.substring(refPosition, refPosition+length);
    			queAln += querySeq.substring(queryPostion, queryPostion+length);
    			refPosition += length;
    			queryPostion += length;
    		}else if( matcher.group(2).compareTo("I") == 0 ) {
    			int length = Integer.parseInt(matcher.group(1));
    			refAln += StringUtils.repeat("-", length);
    			queAln += querySeq.substring(queryPostion, queryPostion+length);
    			queryPostion += length;
    		}else if( matcher.group(2).compareTo("D") == 0 || matcher.group(2).compareTo("N") == 0 ) {
    			int length = Integer.parseInt(matcher.group(1));
    			refAln += referenceSeq.substring(refPosition, refPosition+length);
    			queAln += StringUtils.repeat("-", length);
    			refPosition += length;
    		}else {
    			System.out.println("some thing unknow happened");
    		}
    	}
    	return new Alignment(refAln, queAln);
	}
	public static void main( String argv[] ) {
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/home/bs674/Desktop/OverLapCnsWithGenomeReplictionRegion/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
		SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/home/bs674/Desktop/OverLapCnsWithGenomeReplictionRegion/5.bam"));
		SAMRecordIterator it = reader.iterator();
		while (it.hasNext()) {
            final SAMRecord samRecord = it.next();
            System.out.println(samRecord.getSAMString());
            Alignment alignment = (new CigarToAlignment()).getAlignment(maizeGenome, samRecord);
            System.out.println( alignment.getReferenceAlignment() );
            System.out.println( alignment.getQueryAlignment() );
		}
	}
}
