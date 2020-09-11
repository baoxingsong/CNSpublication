package me.songbx.mains;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CountCNSLengthFromControspondingGenome {
	public static void main(String[] args) {
		new CountCNSLengthFromControspondingGenome("/media/bs674/1twdblack/allCNSbamFiles/sorghum.bam");
		new CountCNSLengthFromControspondingGenome("/media/bs674/1twdblack/allCNSbamFiles/sugarcane_tareploid.bam");
		new CountCNSLengthFromControspondingGenome("/media/bs674/1twdblack/allCNSbamFiles/Hyparrhenia_diplandra.bam");
		new CountCNSLengthFromControspondingGenome("/media/bs674/1twdblack/allCNSbamFiles/Miscanthus_sinensis.bam");
		new CountCNSLengthFromControspondingGenome("/media/bs674/1twdblack/allCNSbamFiles/Chrysopogonserrulatus.bam");
		
	}
	
	
	public CountCNSLengthFromControspondingGenome( String bamFile ) {
		
		try {
			HashSet<String> validatedMaizeChrs = new HashSet<String>();
			validatedMaizeChrs.add("1");
			validatedMaizeChrs.add("2");
			validatedMaizeChrs.add("3");
			validatedMaizeChrs.add("4");
			validatedMaizeChrs.add("5");
			validatedMaizeChrs.add("6");
			validatedMaizeChrs.add("7");
			validatedMaizeChrs.add("8");
			validatedMaizeChrs.add("9");
			validatedMaizeChrs.add("10");
			
			HashMap<String, HashSet<Integer>> positions = new HashMap<String, HashSet<Integer>>();
			
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
			SAMRecordIterator it = reader0.iterator();
			
			Pattern p = Pattern.compile("^(\\d+)H");
			
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            Matcher m=p.matcher(samRecord.getCigarString());
	            if ( validatedMaizeChrs.contains(samRecord.getContig()) && m.find() ) {
		            int queryStart = Integer.parseInt(m.group(1));
	            	int queryEnd = queryStart + samRecord.getReadString().length()-1;
		            
		            if ( samRecord.getReadNegativeStrandFlag() ) {
	    				for ( int position = queryEnd; position >= queryStart; ) {
	        				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
	        				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
	        				while (matcher.find()) {
	        					int length = Integer.parseInt( matcher.group(1));
	        					char cigerLetter =  matcher.group(2).charAt(0);
	        					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
	        						for( int il=0; il<length; ++il ) {
	        							if ( ! positions.containsKey(samRecord.getReadName()) ) {
	        								positions.put(samRecord.getReadName(), new HashSet<Integer>());
	        							}
	        							positions.get(samRecord.getReadName()).add(position);
	        							--position;
	        						}
	        					} else if ( cigerLetter=='I' ) {
	        						position-=length;
	        					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
	        						
	        					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
	        					
	        					}else {
	        						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
	        					}
	        				}
	        			}
	    			} else {
	    				for ( int position = queryStart; position<=queryEnd; ) {
	        				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
	        				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
	        				while (matcher.find()) {
	        					int length = Integer.parseInt( matcher.group(1));
	        					char cigerLetter =  matcher.group(2).charAt(0);
	        					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
	        						for( int il=0; il<length; ++il ) {
	        							if ( ! positions.containsKey(samRecord.getReadName()) ) {
	        								positions.put(samRecord.getReadName(), new HashSet<Integer>());
	        							}
	        							positions.get(samRecord.getReadName()).add(position);
	        							++position;
	        						}
	        					} else if ( cigerLetter=='I' ) {
	        						position+=length;
	        					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
	        						
	        					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
	        					
	        					}else {
	        						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
	        					}
	        				}
	        			}
	    			}	            
	            }
			}
			reader0.close();
			int size = 0;
			for ( String queryContig : positions.keySet() ) {
				size += positions.get(queryContig).size();
			}
			System.out.println(size);
			
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
	
}
//
//55230761
//256214345
//170727578
//114257830
//81062744

