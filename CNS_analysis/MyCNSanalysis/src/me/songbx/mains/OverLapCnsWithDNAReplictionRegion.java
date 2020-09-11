package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.service.IfIntron;

public class OverLapCnsWithDNAReplictionRegion  {
	class BedGraph implements Comparable<BedGraph>{
		private int start;
		private int end;
		private int tag; // 1: in genetic region   2: overlapped with CNS group 1,    3: overlapped with CNS group 2 
		private double score;
		private int phase;
		public BedGraph(int start, int end, double score, int phase) {
			super();
			this.start = start;
			this.end = end;
			this.score = score;
			this.tag = 0;
			this.phase = phase;
		}
		public int getPhase() {
			return phase;
		}
		public void setPhase(int phase) {
			this.phase = phase;
		}
		public int getStart() {
			return start;
		}
		public void setStart(int start) {
			this.start = start;
		}
		public int getEnd() {
			return end;
		}
		public void setEnd(int end) {
			this.end = end;
		}
		public int getTag() {
			return tag;
		}
		public void setTag(int tag) {
			this.tag = tag;
		}
		public double getScore() {
			return score;
		}
		public void setScore(double score) {
			this.score = score;
		}
		@Override
		public int compareTo(BedGraph arg0) {
			return this.getStart() - arg0.getStart();
		}
	}
	
	public HashMap<String, ArrayList<BedGraph>> readBedGraphFile( String bedGraphFile){
		HashMap<String, ArrayList<BedGraph>> bedGraphs = new HashMap<String, ArrayList<BedGraph>>();
		try {
        	File file = new File(bedGraphFile);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = (tempString.trim()).split("\\s+");
				if( !bedGraphs.containsKey( arrOfStr[0]) ) {
					bedGraphs.put(arrOfStr[0], new ArrayList<BedGraph>());
				}
				BedGraph bedGraph = new BedGraph(Integer.parseInt(arrOfStr[1]), Integer.parseInt(arrOfStr[2]), Double.parseDouble(arrOfStr[3]), Integer.parseInt(arrOfStr[7]));
				bedGraphs.get(arrOfStr[0]).add(bedGraph);
			}
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
		for( String chr : bedGraphs.keySet() ) {
			Collections.sort(bedGraphs.get(chr));
		}
		return bedGraphs;
	}
	
	class BedRecord implements Comparable<BedRecord>{
		private int start;
		private int end;
		private int tag;
		public BedRecord(int start, int end) {
			super();
			this.start = start;
			this.end = end;
			this.tag=0;
		}
		public int getTag() {
			return tag;
		}
		public void setTag(int tag) {
			this.tag = tag;
		}
		public int getStart() {
			return start;
		}
		public void setStart(int start) {
			this.start = start;
		}
		public int getEnd() {
			return end;
		}
		public void setEnd(int end) {
			this.end = end;
		}
		@Override
		public int compareTo(BedRecord arg0) {
			return this.getStart() - arg0.getStart();
		}
	}
	
	public void readBamFile(String bamFile, HashMap<String, ArrayList<BedRecord>> bedRecords) {
		try {
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
			SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            int start = samRecord.getAlignmentStart();
	            int end = samRecord.getAlignmentEnd();
	            String chr = samRecord.getContig();
	            if( bedRecords.containsKey(chr) ) {
	            	BedRecord bedRecord = new BedRecord(start, end);
	            	bedRecords.get(chr).add(bedRecord);
	            }
			}
			reader0.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public OverLapCnsWithDNAReplictionRegion( String bedGraphFile, String outputFile ){
		
		HashMap<String, ArrayList<BedGraph>> bedGraphs =  readBedGraphFile(  bedGraphFile);
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
		HashMap<String, ArrayList<BedRecord>> cnsRecords1 = new HashMap<String, ArrayList<BedRecord>> ();
		HashMap<String, ArrayList<BedRecord>> cnsRecords2 = new HashMap<String, ArrayList<BedRecord>> ();
		HashMap<String, ArrayList<BedRecord>> cnsRecords3 = new HashMap<String, ArrayList<BedRecord>> ();
		
		for( int i=1; i<=10; ++i ) {
			String chr = Integer.toString(i);
			cnsRecords1.put(chr, new ArrayList<BedRecord>());
			cnsRecords2.put(chr, new ArrayList<BedRecord>());
			cnsRecords3.put(chr, new ArrayList<BedRecord>());
		}
		
		readBamFile("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam", cnsRecords1);		
		readBamFile("/media/bs674/1_8t/AndCns/CNSprofiles/loop.sam", cnsRecords2);
		readBamFile("/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.sam", cnsRecords3);
		
		HashMap<Integer, Integer> summary1 = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> summary2 = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> summary3 = new HashMap<Integer, Integer>();
		
		try {
			PrintWriter outPut = new PrintWriter(outputFile);
			for( int j=1; j<=10; ++j ) {
				String chr = Integer.toString(j);
				
				Collections.sort(cnsRecords1.get(chr));
				Collections.sort(cnsRecords2.get(chr));
				Collections.sort(cnsRecords3.get(chr));
				
				for( BedGraph bedGraph : bedGraphs.get(chr) ) {
					for( int i=bedGraph.getStart(); i<bedGraph.getEnd(); ++i ) {
						int elements = ifIntron.getElement(chr, i+1);
						if( elements != 0) {
							bedGraph.setTag(  (bedGraph.getTag() | 1) );
							break;
						}
					}
					for( BedRecord cnsRecord : cnsRecords1.get(chr) ) {
						if ( (bedGraph.getStart() <= cnsRecord.getStart() && cnsRecord.getStart()<= bedGraph.getEnd()) ||
								(bedGraph.getStart() <= cnsRecord.getEnd() && cnsRecord.getEnd()<= bedGraph.getEnd()) ||
										(cnsRecord.getStart() <= bedGraph.getStart() && bedGraph.getStart()<= cnsRecord.getEnd()) ||
										(cnsRecord.getStart() <= bedGraph.getEnd() && bedGraph.getEnd()<= cnsRecord.getEnd())
								){
							bedGraph.setTag(  (bedGraph.getTag() | 2) );
							cnsRecord.setTag(cnsRecord.getTag() | bedGraph.getPhase());
						}
					}
					for( BedRecord cnsRecord : cnsRecords2.get(chr) ) {
						if ( (bedGraph.getStart() <= cnsRecord.getStart() && cnsRecord.getStart()<= bedGraph.getEnd()) ||
								(bedGraph.getStart() <= cnsRecord.getEnd() && cnsRecord.getEnd()<= bedGraph.getEnd()) ||
								(cnsRecord.getStart() <= bedGraph.getStart() && bedGraph.getStart()<= cnsRecord.getEnd()) ||
								(cnsRecord.getStart() <= bedGraph.getEnd() && bedGraph.getEnd()<= cnsRecord.getEnd())
							
								){
							bedGraph.setTag(  (bedGraph.getTag() | 4) );
							cnsRecord.setTag(cnsRecord.getTag() | bedGraph.getPhase());
						}
					}
					for( BedRecord cnsRecord : cnsRecords3.get(chr) ) {
						if ( (bedGraph.getStart() <= cnsRecord.getStart() && cnsRecord.getStart()<= bedGraph.getEnd()) ||
								(bedGraph.getStart() <= cnsRecord.getEnd() && cnsRecord.getEnd()<= bedGraph.getEnd()) ||
								(cnsRecord.getStart() <= bedGraph.getStart() && bedGraph.getStart()<= cnsRecord.getEnd()) ||
								(cnsRecord.getStart() <= bedGraph.getEnd() && bedGraph.getEnd()<= cnsRecord.getEnd())
								){
							bedGraph.setTag(  (bedGraph.getTag() | 8) );
							cnsRecord.setTag(cnsRecord.getTag() | bedGraph.getPhase());
						}
					}
					outPut.println(chr + "\t" + bedGraph.getStart()+"\t"+bedGraph.getEnd()+"\t"+bedGraph.getScore()+"\t"+bedGraph.getTag());
				}
				
				for( BedRecord cnsRecord : cnsRecords1.get(chr) ) {
					if( ! summary1.containsKey(cnsRecord.getTag()) ) {
						summary1.put(cnsRecord.getTag(), 0);
					}
					summary1.put(cnsRecord.getTag(), summary1.get(cnsRecord.getTag())+1);
				}
				
				for( BedRecord cnsRecord : cnsRecords2.get(chr) ) {
					if( ! summary2.containsKey(cnsRecord.getTag()) ) {
						summary2.put(cnsRecord.getTag(), 0);
					}
					summary2.put(cnsRecord.getTag(), summary2.get(cnsRecord.getTag())+1);
				}
				for( BedRecord cnsRecord : cnsRecords3.get(chr) ) {
					if( ! summary3.containsKey(cnsRecord.getTag()) ) {
						summary3.put(cnsRecord.getTag(), 0);
					}
					summary3.put(cnsRecord.getTag(), summary3.get(cnsRecord.getTag())+1);
				}
			}
			outPut.println("#summary1");
			for( int key : summary1.keySet() ) {
				outPut.println("#"+key + "\t" + summary1.get(key));
			}
			outPut.println("#summary2");
			for( int key : summary2.keySet() ) {
				outPut.println("#"+key + "\t" + summary2.get(key));
			}
			outPut.println("#summary3");
			for( int key : summary3.keySet() ) {
				outPut.println("#"+key + "\t" + summary3.get(key));
			}
			outPut.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		new OverLapCnsWithDNAReplictionRegion("/media/bs674/1_8t/AndCns/overlapWithGenomeReplication/Repli_S.table", "/media/bs674/1_8t/AndCns/overlapWithGenomeReplication/Repli_S_tag.bedgraph");
	}
}
