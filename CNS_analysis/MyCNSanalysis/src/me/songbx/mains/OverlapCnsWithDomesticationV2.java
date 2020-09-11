package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.model.BedRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class OverlapCnsWithDomesticationV2 {
	public static void main1(String[] args) {
		try {
			HashMap<String, ArrayList<BedRecord>> transcripts = new HashMap<String, ArrayList<BedRecord>>();
			{
				File file = new File("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = reader.readLine(); // escape the first line
				Pattern p = Pattern.compile("^(\\S+)\t\\S+\tmRNA\t(\\d+)\t(\\d+).*biotype=protein_coding");
				while ((tempString = reader.readLine()) != null) {
					Matcher m=p.matcher(tempString);
					if(m.find()){
						String chr = m.group(1);
						int start = Integer.parseInt(m.group(2));
						int end = Integer.parseInt(m.group(3));
			            if( transcripts.containsKey(chr) ) {
			            	
			            }else {
			            	transcripts.put(chr, new ArrayList<BedRecord>());
			            }
			            transcripts.get(chr).add(new BedRecord(start, end));
					}
				}
				reader.close();
			}
			System.out.println("gff file reading done");
			
			
			HashMap<String, ArrayList<BedRecord>> cnss = new HashMap<String, ArrayList<BedRecord>>();
			
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
			SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            String chr = samRecord.getContig();
	            int start = samRecord.getAlignmentStart();
	            int end = samRecord.getAlignmentEnd();
	            if( cnss.containsKey(chr) ) {
	            	
	            }else {
	            	cnss.put(chr, new ArrayList<BedRecord>());
	            }
	            cnss.get(chr).add(new BedRecord(start, end));
			}
			reader0.close();
			
			{
				PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/overLapWithDomestication/localAdaptation/RAISD_with_distance");
				File file = new File("/media/bs674/1_8t/AndCns/overLapWithDomestication/localAdaptation/FullRAISD_outliers.csv");
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = reader.readLine(); // escape the first line
				Pattern p = Pattern.compile("^(\\S+)\t(\\d+)\t(\\S+)\t");
				while ((tempString = reader.readLine()) != null) {
					Matcher m=p.matcher(tempString);
					if(m.find()){
						String chr = m.group(1);
						int position = Integer.parseInt(m.group(2));
						boolean withinTranscript = false;
						
						for( BedRecord transcript : transcripts.get(chr) ) {
							if( position>=transcript.getStart() && position<=transcript.getEnd()  ) {
								withinTranscript = true;
							}
						}
						
						if( !withinTranscript ) {
							int miniDistance = 1000000000;
							for( BedRecord bedRecord : cnss.get(chr) ) {
								if( position>=bedRecord.getStart() && position<=bedRecord.getEnd()  ) {
									miniDistance = 0;
								}else {
									miniDistance = miniDistance < Math.abs(bedRecord.getStart()-position) ? miniDistance : Math.abs(bedRecord.getStart()-position);
									miniDistance = miniDistance < Math.abs(bedRecord.getEnd()-position) ? miniDistance : Math.abs(bedRecord.getEnd()-position);
								}
							}
							outPut.println(tempString + "\t" + miniDistance);
						}else {
							outPut.println(tempString + "\t-10");
						}
					}else {
						System.err.println(tempString);
					}
				}
				outPut.close();
				reader.close();
			}
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static void main(String[] args) {
		try {
			HashMap<String, ArrayList<BedRecord>> transcripts = new HashMap<String, ArrayList<BedRecord>>();
			{
				File file = new File("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = reader.readLine(); // escape the first line
				Pattern p = Pattern.compile("^(\\S+)\t\\S+\tmRNA\t(\\d+)\t(\\d+).*biotype=protein_coding");
				while ((tempString = reader.readLine()) != null) {
					Matcher m=p.matcher(tempString);
					if(m.find()){
						String chr = m.group(1);
						int start = Integer.parseInt(m.group(2));
						int end = Integer.parseInt(m.group(3));
			            if( transcripts.containsKey(chr) ) {
			            	
			            }else {
			            	transcripts.put(chr, new ArrayList<BedRecord>());
			            }
			            transcripts.get(chr).add(new BedRecord(start, end));
					}
				}
				reader.close();
			}
			System.out.println("gff file reading done");
			
			{
				PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/overLapWithDomestication/localAdaptation/RAISD_with_distance_to_transcript");
				File file = new File("/media/bs674/1_8t/AndCns/overLapWithDomestication/localAdaptation/FullRAISD_outliers.csv");
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = reader.readLine(); // escape the first line
				Pattern p = Pattern.compile("^(\\S+)\t(\\d+)\t(\\S+)\t");
				while ((tempString = reader.readLine()) != null) {
					Matcher m=p.matcher(tempString);
					if(m.find()){
						String chr = m.group(1);
						int position = Integer.parseInt(m.group(2));
					
						int miniDistance = 1000000000;
						for( BedRecord bedRecord : transcripts.get(chr) ) {
							if( position>=bedRecord.getStart() && position<=bedRecord.getEnd()  ) {
								miniDistance = 0;
							}else {
								miniDistance = miniDistance < Math.abs(bedRecord.getStart()-position) ? miniDistance : Math.abs(bedRecord.getStart()-position);
								miniDistance = miniDistance < Math.abs(bedRecord.getEnd()-position) ? miniDistance : Math.abs(bedRecord.getEnd()-position);
							}
						}
						outPut.println(tempString + "\t" + miniDistance);
					
					}else {
						System.err.println(tempString);
					}
				}
				outPut.close();
				reader.close();
			}
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}
