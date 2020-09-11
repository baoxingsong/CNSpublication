package me.songbx.service;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.impl.AnnotationReadImpl;
import me.songbx.model.Cds;
import me.songbx.model.GeneSimple;
import me.songbx.model.NearByGene;
import me.songbx.model.Strand;
import me.songbx.model.Transcript;

public class IfIntron {
	private HashMap<String, ArrayList<Transcript>> transcriptHashMap = new HashMap<String, ArrayList<Transcript>>();
	private HashMap<String, ArrayList<GeneSimple>> geneHashMap = new HashMap<String, ArrayList<GeneSimple>>();
	public HashMap<String, ArrayList<Transcript>> getTranscriptHashMap() {
		return transcriptHashMap;
	}
	public void setTranscriptHashMap(HashMap<String, ArrayList<Transcript>> transcriptHashMap) {
		this.transcriptHashMap = transcriptHashMap;
	}
	public HashMap<String, ArrayList<GeneSimple>> getGeneHashMap() {
		return geneHashMap;
	}

	public void setGeneHashMap(HashMap<String, ArrayList<GeneSimple>> geneHashMap) {
		this.geneHashMap = geneHashMap;
	}



	public IfIntron( String fileLocation ) {
		HashMap<String, HashSet<Transcript>> transcriptHashSet = new AnnotationReadImpl(fileLocation).getTranscriptHashSet();
		for ( String chr : transcriptHashSet.keySet() ) {
			transcriptHashMap.put(chr, new ArrayList<Transcript>());
			for( Transcript transcript : transcriptHashSet.get(chr) ){
				transcriptHashMap.get(chr).add(transcript);
			}
			Collections.sort(transcriptHashMap.get(chr));
		}
		File file = new File(fileLocation);
		BufferedReader reader = null;
        try {
        	reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            Pattern p = Pattern.compile("^(\\S+)\\s+\\S+\\s+mRNA\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\S+)\\s+.*;Parent=gene:(\\w+);.*protein_coding");
			while ((tempString = reader.readLine()) != null) {
            	Matcher m=p.matcher(tempString);
				if(tempString.startsWith("#")){
					
				}else{
					Matcher mf = null;
					if(m.find()){
						mf = m;
					}
	                if(mf != null){
	                	int start = Integer.parseInt(mf.group(2));
	                	int end = Integer.parseInt(mf.group(3));
	                	if(start>end){
	                		int temp=start;
	                		start = end;
	                		end = temp;
	                	}
	                	Strand strand;
	                	if("-".equals(mf.group(4))){
	                		strand = Strand.NEGTIVE;
	                	}else{
	                		strand = Strand.POSITIVE;
	                	}
	                	String chrId = mf.group(1);
	                	String geneId = mf.group(5);
	                	GeneSimple geneSimple = new GeneSimple(strand, chrId, start, end, geneId);
                		if( ! geneHashMap.containsKey(chrId) ) {
                			geneHashMap.put(chrId, new ArrayList<GeneSimple>());
                		}
                		geneHashMap.get(chrId).add(geneSimple);
	                }
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {
                	e1.getStackTrace();
                }
            }
        }
        for ( String chr : geneHashMap.keySet() ) {
			Collections.sort(geneHashMap.get(chr));
		}
	}
	
	
	
	// the input a coordinate
	// the input is a int, 
	// 0 is intergenic
	// 1 is intron
	// 2 is CDS
	// 3 UTR
	public int getElement ( String chr, int position ) {  // find the overlaped gene using a half interval search approach
		int element = 0; //intergenic
		if( this.geneHashMap.containsKey(chr) ) {
			int start = 0 ;
			int end = geneHashMap.get(chr).size()-1;
			int lastStart = start;
			if( geneHashMap.get(chr).get(start).getStart() > position ){
				return element;
			}
			if(geneHashMap.get(chr).get(end).getEnd() < position){
				return element;
			}
//			System.out.println("line 101");
			if( end > start ) {
//				System.out.println("line 103:" + start + "\t" + end);
				while(!((geneHashMap.get(chr).get(start).getStart() <= position) && (geneHashMap.get(chr).get(start+1).getEnd( )>= position))){
//					System.out.println("line 105:" + start + "\t" + end);
					if((geneHashMap.get(chr).get(start).getEnd() < position)){
//						System.out.println("line 107:" + start + "\t" + end);
						lastStart = start;
						if(1 == (end - start)){
							start = end;
						}else{
							start = (start+end)/2;
						}
					}else{
//						System.out.println("line 115:" + start + "\t" + end);
						end = start;
						start = lastStart;
						if( start == end ) {
							break;
						}
					}
				}
			}
			if( geneHashMap.get(chr).get(start).getStart()<=position && position<=geneHashMap.get(chr).get(start).getEnd()) {
				element = 3;
			} else if ( geneHashMap.get(chr).get(start+1).getStart()<=position && position<=geneHashMap.get(chr).get(start+1).getEnd() ) {
				element = 3;
			}
//			System.out.println("line 123");
			if( 3==element ) {
				if( this.transcriptHashMap.containsKey(chr) ) {
					start = 0 ;
					end = transcriptHashMap.get(chr).size()-1;
					lastStart = start;
					if( transcriptHashMap.get(chr).get(start).getStart() > position ){
						return element;
					}
					if(transcriptHashMap.get(chr).get(end).getEnd() < position){
						return element;
					}
//					System.out.println("line 135");
					if( end > start ) {
//						System.out.println("line 137:" + start + "\t" + end);
						while(!((transcriptHashMap.get(chr).get(start).getStart() <= position) && (transcriptHashMap.get(chr).get(start+1).getEnd( )>= position))){
							if((transcriptHashMap.get(chr).get(start).getEnd() < position)){
								lastStart = start;
								if(1 == (end - start)){
									start = end;
								}else{
									start = (start+end)/2;
								}
							}else{
								end = start;
								start = lastStart;
								if( start == end ) {
									break;
								}
							}
						}
					}
					if( transcriptHashMap.get(chr).get(start).getStart()<=position && position<=transcriptHashMap.get(chr).get(start).getEnd()) {
						element = 1;
						for ( Cds cds : transcriptHashMap.get(chr).get(start).getCdsHashSet()  ) {
							if ( cds.getStart()<= position && position <=cds.getEnd()  ) {
								element = 2;
							}
						}
					}
					if ( transcriptHashMap.get(chr).size()>1 && transcriptHashMap.get(chr).get(start+1).getStart()<=position && position<=transcriptHashMap.get(chr).get(start+1).getEnd() ) {
						element = 1;
						for ( Cds cds : transcriptHashMap.get(chr).get(start+1).getCdsHashSet()  ) {
							if ( cds.getStart()<= position && position <=cds.getEnd()  ) {
								element = 2;
							}
						}
					}			
				}
			}
		}
		return element;
	}
	
	public NearByGene getMinidistanceWithTranscript ( String chr, int position1, int position2 ) {
		int miniDistance = Integer.MAX_VALUE;
		String gene="";
		if( this.geneHashMap.containsKey(chr) ) {
			for( GeneSimple transcript : geneHashMap.get(chr) ) {
				int distance1 = Integer.MAX_VALUE;
				int distance2 = Integer.MAX_VALUE;
				if( Strand.NEGTIVE == transcript.getStrand() ) {
					if( position1 < transcript.getStart() ) {
						distance1 = (transcript.getStart() - position1); //positive 3'
					}else if ( position1 > transcript.getEnd() ) {
						distance1 = -(position1 - transcript.getEnd()); //negative 5'
					}
					
					if( position2 < transcript.getStart() ) {
						distance2 = (transcript.getStart() - position2);
					}else if ( position2 > transcript.getEnd() ) {
						distance2 = -(position2 - transcript.getEnd());
					}
				}else {
					if( position1 < transcript.getStart() ) {
						distance1 = -(transcript.getStart() - position1);//negative 5'
					}else if ( position1 > transcript.getEnd() ) {
						distance1 = (position1 - transcript.getEnd());
					}
					if( position2 < transcript.getStart() ) {
						distance2 = -(transcript.getStart() - position2);
					}else if ( position2 > transcript.getEnd() ) {
						distance2 = (position2 - transcript.getEnd());
					}
				}
				if( Math.abs(miniDistance) > Math.abs(distance1) ) {
					miniDistance = distance1;
					gene=transcript.getName();
				}
				if( Math.abs(miniDistance) > Math.abs(distance2) ) {
					miniDistance = distance2;
					gene=transcript.getName();
				}
			}
		}
		return (new NearByGene(miniDistance, gene));
	}
	
	public NearByGene getMinidistanceWithTranscript ( String chr, int position ) {
		int miniDistance = Integer.MAX_VALUE;
		String gene="";
		if( this.geneHashMap.containsKey(chr) ) {
			for( GeneSimple transcript : geneHashMap.get(chr) ) {
				int distance = Integer.MAX_VALUE;
				if( Strand.NEGTIVE == transcript.getStrand() ) {
					if( position < transcript.getStart() ) {
						distance = (transcript.getStart() - position); //positive 3'
					}else if ( position > transcript.getEnd() ) {
						distance = -(position - transcript.getEnd()); //negative 5'
					}
				}else {
					if( position < transcript.getStart() ) {
						distance = -(transcript.getStart() - position);//negative 5'
					}else if ( position > transcript.getEnd() ) {
						distance = (position - transcript.getEnd());
					}
				}
				if( Math.abs(miniDistance) > Math.abs(distance) ) {
					miniDistance = distance;
					gene=transcript.getName();
				}
			}
		}
		return (new NearByGene(miniDistance, gene));
	}
	
	public NearByGene getMinidistanceWithCds ( String chr, int position ) {
		int miniDistance = Integer.MAX_VALUE;
		String gene="";
		if( this.transcriptHashMap.containsKey(chr) ) {
			for( Transcript transcript : transcriptHashMap.get(chr) ) {
				int distance = Integer.MAX_VALUE;
				if( Strand.NEGTIVE == transcript.getStrand() ) {
					if( position < transcript.getStart() ) {
						distance = (transcript.getStart() - position); //positive 3'
					}else if ( position > transcript.getEnd() ) {
						distance = -(position - transcript.getEnd());  //negative 5'
					}
				}else {
					if( position < transcript.getStart() ) {
						distance = -(transcript.getStart() - position); //negative 5'
					}else if ( position > transcript.getEnd() ) {
						distance = (position - transcript.getEnd()); //positive 3'
					}
				}
				if( Math.abs(miniDistance) > Math.abs(distance) ) {
					miniDistance = distance;
					gene=transcript.getName();
				}
			}
		}
		return (new NearByGene(miniDistance, gene));
	}
	
	public NearByGene getMinidistanceWithCds ( String chr, int position1, int position2 ) {
		int miniDistance = Integer.MAX_VALUE;
		String gene="";
		if( this.transcriptHashMap.containsKey(chr) ) {
			for( Transcript transcript : transcriptHashMap.get(chr) ) {
				int distance1 = Integer.MAX_VALUE;
				int distance2 = Integer.MAX_VALUE;
				if( Strand.NEGTIVE == transcript.getStrand() ) {
					if( position1 < transcript.getStart() ) {
						distance1 = (transcript.getStart() - position1); //positive 3'
					}else if ( position1 > transcript.getEnd() ) {
						distance1 = -(position1 - transcript.getEnd());  //negative 5'
					}
					if( position2 < transcript.getStart() ) {
						distance2 = (transcript.getStart() - position2); //positive 3'
					}else if ( position2 > transcript.getEnd() ) {
						distance2 = -(position2 - transcript.getEnd());  //negative 5'
					}
				}else {
					if( position1 < transcript.getStart() ) {
						distance1 = -(transcript.getStart() - position1); //negative 5'
					}else if ( position1 > transcript.getEnd() ) {
						distance1 = (position1 - transcript.getEnd()); //positive 3'
					}
					if( position2 < transcript.getStart() ) {
						distance2 = -(transcript.getStart() - position2); //negative 5'
					}else if ( position2 > transcript.getEnd() ) {
						distance2 = (position2 - transcript.getEnd()); //positive 3'
					}
				}
				if( Math.abs(miniDistance) > Math.abs(distance1) ) {
					miniDistance = distance1;
					gene=transcript.getName();
				}
				if( Math.abs(miniDistance) > Math.abs(distance2) ) {
					miniDistance = distance2;
					gene=transcript.getName();
				}
			}
		}
		return (new NearByGene(miniDistance, gene));
	}
	
	public static void main( String[] args ) {
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/sorghum/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.42.gff3");
		System.out.println("reading done");
		System.out.println(ifIntron.getElement("6", 16538)); //1
		System.out.println(ifIntron.getElement("6", 12300)); //2 there are two genes located at opposting strand, very complex
//		System.out.println(ifIntron.getElement("1", 71982186)); //2
//		System.out.println(ifIntron.getElement("1", 71982187)); //1
//		System.out.println(ifIntron.getElement("1", 71981858)); //0
		
		HashMap<String, ArrayList<Transcript>> transcriptHashMap = ifIntron.getTranscriptHashMap();
		HashMap<String, ArrayList<GeneSimple>> geneHashMap = ifIntron.getGeneHashMap();
		HashMap<String, HashSet<Integer>> allGeneticsBps = new HashMap<String, HashSet<Integer>>();
		HashMap<String, HashSet<Integer>> allCdsBps = new HashMap<String, HashSet<Integer>>();
		int allCdsLength=0;
		int allGeneticLength=0;
		for( int i=1; i<=10; i++ ) {
			String chr = Integer.toString(i);
			allGeneticsBps.put(chr, new HashSet<Integer>());
			for( GeneSimple geneSimple : geneHashMap.get(chr) ) {
				for( int position=geneSimple.getStart(); position<=geneSimple.getEnd(); ++position ) {
					allGeneticsBps.get(chr).add(position);
				}
			}
			allCdsBps.put(chr, new HashSet<Integer>());
			for( Transcript transcript : transcriptHashMap.get(chr) ) {
				for( Cds cds :  transcript.getCdsHashSet()) {
					for( int position=cds.getStart(); position<=cds.getEnd(); ++position ) {
						allCdsBps.get(chr).add(position);
					}
				}
			}
			allCdsLength += allCdsBps.get(chr).size();
			allGeneticLength += allGeneticsBps.get(chr).size();
		}
		System.out.println("allCdsLength:" + allCdsLength);
		System.out.println("allGeneticLength:" + allGeneticLength);
	}
}
