package me.songbx.mains;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.ReadGffForSimpleGene;
import me.songbx.model.GeneSimple;
import me.songbx.model.Strand;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.commons.lang3.StringUtils;

public class CnsCountBpsWithExpressionStablity {
	public static void main(String[] args) {
		HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
		HashMap<String, String> maize2 = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46_cds.fa").getChromoSomeHashMap();
		
		PrintWriter outPut = null;

		int maizeExtendLength = 5000;
		int miniscore= 0;
		int maxscore= 50000;
		
		try {
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/subfunctionlization/sorghum/countofCNSbps"+miniscore+"_"+maxscore+"_"+maizeExtendLength);
			
	        Pattern p = Pattern.compile("^(\\d+)H");
	        Pattern pattern = Pattern.compile("(\\d+)(\\w)");
	        { // all the cis region
	        	
		        for ( String maizeGeneId : maizeGenes.keySet() ) {
		        	HashSet<Integer> allPositions = new HashSet<Integer>();
		        	String maizeChr = maizeGenes.get(maizeGeneId).getChromeSomeName();
					int maizeStart = maizeGenes.get(maizeGeneId).getStart() - maizeExtendLength;
					int maizeEnd = maizeGenes.get(maizeGeneId).getEnd() + maizeExtendLength;
					
					SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
					
					SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
					while (it.hasNext()) {
			            final SAMRecord samRecord = it.next();
			            Matcher m=p.matcher(samRecord.getCigarString());
			            if(m.find() && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore ){
			            	String seqMaize = maize2.get(maizeChr).substring(samRecord.getAlignmentStart()-1, samRecord.getAlignmentEnd());
	            			String seqsorghum = samRecord.getReadString();
	            			String alnMaize = "";
	            			String alnSorghum = "";
	            			int maizePosition = 0;
	            			int sorghumPostion = 0;
	            			Matcher matcher = pattern.matcher(samRecord.getCigarString());
			            	while( matcher.find() ) {
			            		if( matcher.group(2).compareTo("H") == 0 ) {
			            			
			            		}else if( matcher.group(2).compareTo("M") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
			            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
			            			maizePosition += length;
			            			sorghumPostion += length;
			            		}else if( matcher.group(2).compareTo("I") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += StringUtils.repeat("-", length);
			            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
			            			sorghumPostion += length;
			            		}else if( matcher.group(2).compareTo("D") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
			            			alnSorghum += StringUtils.repeat("-", length);
			            			maizePosition += length;
			            		}else {
			            			System.out.println("some thing unknow happened");
			            		}
			            	}
		            		int position = maizeStart;
		            		for( int ip = 0; ip<alnMaize.length(); ++ip ) {
								if( alnMaize.charAt(ip) != 'N' && alnMaize.charAt(ip)==alnSorghum.charAt(ip)  )  {
									allPositions.add(position);
								}
								if( alnMaize.charAt(ip)  != '-') {
									++position;
								}
							}	
	            		}
					}
	        		outPut.println("all\t" + maizeGeneId + "\t" + allPositions.size());
	        		reader.close();
	        	}
	        }
	    
	        { // up stream
		        for ( String maizeGeneId : maizeGenes.keySet() ) {
		        	HashSet<Integer> allPositions = new HashSet<Integer>();
		        	String maizeChr = maizeGenes.get(maizeGeneId).getChromeSomeName();
					int maizeStart = maizeGenes.get(maizeGeneId).getStart() - maizeExtendLength;
					int maizeEnd = maizeGenes.get(maizeGeneId).getStart() - 1;
					
					if ( maizeGenes.get(maizeGeneId).getStrand() == Strand.NEGTIVE ) {
						maizeStart = maizeGenes.get(maizeGeneId).getEnd()+1;
						maizeEnd = maizeGenes.get(maizeGeneId).getEnd()+maizeExtendLength;
					}
					
					SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
					
					SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
					while (it.hasNext()) {
			            final SAMRecord samRecord = it.next();
			            Matcher m=p.matcher(samRecord.getCigarString());
			            if(m.find() && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore ){
			            	String seqMaize = maize2.get(maizeChr).substring(samRecord.getAlignmentStart()-1, samRecord.getAlignmentEnd());
	            			String seqsorghum = samRecord.getReadString();
	            			String alnMaize = "";
	            			String alnSorghum = "";
	            			int maizePosition = 0;
	            			int sorghumPostion = 0;
	            			Matcher matcher = pattern.matcher(samRecord.getCigarString());
			            	while( matcher.find() ) {
			            		if( matcher.group(2).compareTo("H") == 0 ) {
			            			
			            		}else if( matcher.group(2).compareTo("M") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
			            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
			            			maizePosition += length;
			            			sorghumPostion += length;
			            		}else if( matcher.group(2).compareTo("I") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += StringUtils.repeat("-", length);
			            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
			            			sorghumPostion += length;
			            		}else if( matcher.group(2).compareTo("D") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
			            			alnSorghum += StringUtils.repeat("-", length);
			            			maizePosition += length;
			            		}else {
			            			System.out.println("some thing unknow happened");
			            		}
			            	}
		            		int position = maizeStart;
		            		for( int ip = 0; ip<alnMaize.length(); ++ip ) {
								if( alnMaize.charAt(ip) != 'N' && alnMaize.charAt(ip)==alnSorghum.charAt(ip)  )  {
									allPositions.add(position);
								}
								if( alnMaize.charAt(ip)  != '-') {
									++position;
								}
							}
	            		}
					}
					outPut.println("upstream\t" + maizeGeneId + "\t" + allPositions.size());
					reader.close();
	        	}
	    	}
	        { // inner
	        	for ( String maizeGeneId : maizeGenes.keySet() ) {
	        		if( maizeGenes.get(maizeGeneId).getNumberOfCds() > 1 ) { // inner
	        			HashSet<Integer> allPositions = new HashSet<Integer>();
	        			String maizeChr = maizeGenes.get(maizeGeneId).getChromeSomeName();
			        
	        			int maizeStart = maizeGenes.get(maizeGeneId).getStart();
						int maizeEnd = maizeGenes.get(maizeGeneId).getEnd();
						SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
						
						SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
						while (it.hasNext()) {
				            final SAMRecord samRecord = it.next();
				            Matcher m=p.matcher(samRecord.getCigarString());
				            if(m.find() && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore ){
				            	String seqMaize = maize2.get(maizeChr).substring(samRecord.getAlignmentStart()-1, samRecord.getAlignmentEnd());
		            			String seqsorghum = samRecord.getReadString();
		            			String alnMaize = "";
		            			String alnSorghum = "";
		            			int maizePosition = 0;
		            			int sorghumPostion = 0;
		            			Matcher matcher = pattern.matcher(samRecord.getCigarString());
				            	while( matcher.find() ) {
				            		if( matcher.group(2).compareTo("H") == 0 ) {
				            		}else if( matcher.group(2).compareTo("M") == 0 ) {
				            			int length = Integer.parseInt(matcher.group(1));
				            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
				            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
				            			maizePosition += length;
				            			sorghumPostion += length;
				            		}else if( matcher.group(2).compareTo("I") == 0 ) {
				            			int length = Integer.parseInt(matcher.group(1));
				            			alnMaize += StringUtils.repeat("-", length);
				            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
				            			sorghumPostion += length;
				            		}else if( matcher.group(2).compareTo("D") == 0 ) {
				            			int length = Integer.parseInt(matcher.group(1));
				            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
				            			alnSorghum += StringUtils.repeat("-", length);
				            			maizePosition += length;
				            		}else {
				            			System.out.println("some thing unknow happened");
				            		}
				            	}
			            		int position = maizeStart;
			            		for( int ip = 0; ip<alnMaize.length(); ++ip ) {
									if( alnMaize.charAt(ip) != 'N' && alnMaize.charAt(ip)==alnSorghum.charAt(ip)  )  {
										allPositions.add(position);
									}
									if( alnMaize.charAt(ip)  != '-') {
										++position;
									}
								}
		            		}
						}
						outPut.println("intron\t" + maizeGeneId + "\t" + allPositions.size());
						reader.close();		
	        		}
	        	}
	        }
								
			{ // down stream
				for ( String maizeGeneId : maizeGenes.keySet() ) {
		        	HashSet<Integer> allPositions = new HashSet<Integer>();
		        	String maizeChr = maizeGenes.get(maizeGeneId).getChromeSomeName();
		        	
					int maizeStart = maizeGenes.get(maizeGeneId).getEnd()+1;
					int maizeEnd = maizeGenes.get(maizeGeneId).getEnd()+maizeExtendLength;
					
					if ( maizeGenes.get(maizeGeneId).getStrand() == Strand.NEGTIVE ) {
						maizeStart = maizeGenes.get(maizeGeneId).getStart() - maizeExtendLength;
						maizeEnd = maizeGenes.get(maizeGeneId).getStart()-1;
					}
					
					SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
					
					SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
					while (it.hasNext()) {
			            final SAMRecord samRecord = it.next();
			            Matcher m=p.matcher(samRecord.getCigarString());
			            if(m.find() && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore ){
			            	String seqMaize = maize2.get(maizeChr).substring(samRecord.getAlignmentStart()-1, samRecord.getAlignmentEnd());
	            			String seqsorghum = samRecord.getReadString();
	            			String alnMaize = "";
	            			String alnSorghum = "";
	            			int maizePosition = 0;
	            			int sorghumPostion = 0;
	            			Matcher matcher = pattern.matcher(samRecord.getCigarString());
			            	while( matcher.find() ) {
			            		if( matcher.group(2).compareTo("H") == 0 ) {
			            			
			            		}else if( matcher.group(2).compareTo("M") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
			            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
			            			maizePosition += length;
			            			sorghumPostion += length;
			            		}else if( matcher.group(2).compareTo("I") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += StringUtils.repeat("-", length);
			            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
			            			sorghumPostion += length;
			            		}else if( matcher.group(2).compareTo("D") == 0 ) {
			            			int length = Integer.parseInt(matcher.group(1));
			            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
			            			alnSorghum += StringUtils.repeat("-", length);
			            			maizePosition += length;
			            		}else {
			            			System.out.println("some thing unknow happened");
			            		}
			            	}
		            		int position = maizeStart;
		            		for( int ip = 0; ip<alnMaize.length(); ++ip ) {
								if( alnMaize.charAt(ip) != 'N' && alnMaize.charAt(ip)==alnSorghum.charAt(ip)  )  {
									allPositions.add(position);
								}
								if( alnSorghum.charAt(ip)  != '-') {
									++position;
								}
							}
	            		}
					}
					outPut.println("downstream\t" + maizeGeneId + "\t" + allPositions.size());
					reader.close();
	        	}
			}
			outPut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
}
