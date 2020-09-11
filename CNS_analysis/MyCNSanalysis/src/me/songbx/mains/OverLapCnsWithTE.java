package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.impl.ChromoSomeReadImpl;

public class OverLapCnsWithTE {
	public static void main(String[] args) {
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		new OverLapCnsWithTE("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.wig", "/media/bs674/1_8t/AndCns/overlapCNSWithDIfferentTEsuperfamily/B73.structuralTEv2.fulllength.2018-09-19.gff3", maizeGenome);
		new OverLapCnsWithTE("/media/bs674/1_8t/AndCns/CNSprofiles/loop.wig", "/media/bs674/1_8t/AndCns/overlapCNSWithDIfferentTEsuperfamily/B73.structuralTEv2.fulllength.2018-09-19.gff3", maizeGenome);
		new OverLapCnsWithTE("/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.wig", "/media/bs674/1_8t/AndCns/overlapCNSWithDIfferentTEsuperfamily/B73.structuralTEv2.fulllength.2018-09-19.gff3", maizeGenome);
		
		
		/*
		 *
teClass enrichmentFold  cnsTeLength     cnsLength       teLength        genomeSize
DHH     1.2651535376222418      1626219.0       2.5708071E7     1.05316009E8    2.106338117E9
DTT     0.38271926610259055     109051.0        2.5708071E7     2.3345755E7     2.106338117E9
RLC     0.07113072628650152     510901.0        2.5708071E7     5.8848867E8     2.106338117E9
DTX     0.03571361079213099     7742.0  2.5708071E7     1.7761433E7     2.106338117E9
RLG     0.03669023034553889     419683.0        2.5708071E7     9.37194078E8    2.106338117E9
DTA     0.1493218870721347      7612.0  2.5708071E7     4176706.0       2.106338117E9
DTC     0.02046999241656658     2482.0  2.5708071E7     9934424.0       2.106338117E9
RST     1.2440999544218858      4314.0  2.5708071E7     284108.0        2.106338117E9
RIL     0.06481463180994428     1715.0  2.5708071E7     2167952.0       2.106338117E9
DTH     0.21559434563359559     55512.0 2.5708071E7     2.1096388E7     2.106338117E9
DTM     0.09679683560962095     669.0   2.5708071E7     566270.0        2.106338117E9
RIT     0.010668984422068085    66.0    2.5708071E7     506850.0        2.106338117E9
RLX     0.2403106736044308      773427.0        2.5708071E7     2.63696803E8    2.106338117E9
teClass enrichmentFold  cnsTeLength     cnsLength       teLength        genomeSize
DHH     0.0     0.0     6964524.0       1.05316009E8    2.106338117E9
DTT     1.1886751861882514      91756.0 6964524.0       2.3345755E7     2.106338117E9
RLC     0.08413854011907047     163718.0        6964524.0       5.8848867E8     2.106338117E9
DTX     0.08347029581815987     4902.0  6964524.0       1.7761433E7     2.106338117E9
RLG     0.042056989103267724    130326.0        6964524.0       9.37194078E8    2.106338117E9
DTA     0.705714675059103       9746.0  6964524.0       4176706.0       2.106338117E9
DTC     0.22543379329345356     7405.0  6964524.0       9934424.0       2.106338117E9
RST     0.0     0.0     6964524.0       284108.0        2.106338117E9
RIL     0.0     0.0     6964524.0       2167952.0       2.106338117E9
DTH     0.5674195999276086      39580.0 6964524.0       2.1096388E7     2.106338117E9
DTM     0.0     0.0     6964524.0       566270.0        2.106338117E9
RIT     0.28880357070821805     484.0   6964524.0       506850.0        2.106338117E9
RLX     0.27377698586370375     238707.0        6964524.0       2.63696803E8    2.106338117E9
teClass enrichmentFold  cnsTeLength     cnsLength       teLength        genomeSize
DHH     0.0     0.0     1.9379352E7     1.05316009E8    2.106338117E9
DTT     0.8390284130705512      180217.0        1.9379352E7     2.3345755E7     2.106338117E9
RLC     0.24965653086616735     1351737.0       1.9379352E7     5.8848867E8     2.106338117E9
DTX     0.12355737204152328     20191.0 1.9379352E7     1.7761433E7     2.106338117E9
RLG     0.2260291512774395      1948970.0       1.9379352E7     9.37194078E8    2.106338117E9
DTA     1.4903809096357639      57272.0 1.9379352E7     4176706.0       2.106338117E9
DTC     1.544578883119123       141177.0        1.9379352E7     9934424.0       2.106338117E9
RST     0.0     0.0     1.9379352E7     284108.0        2.106338117E9
RIL     0.0     0.0     1.9379352E7     2167952.0       2.106338117E9
DTH     0.48905394005601566     94924.0 1.9379352E7     2.1096388E7     2.106338117E9
DTM     0.0     0.0     1.9379352E7     566270.0        2.106338117E9
RIT     0.2792031856540883      1302.0  1.9379352E7     506850.0        2.106338117E9
RLX     0.9065841773291028      2199501.0       1.9379352E7     2.63696803E8    2.106338117E9

		 * */
	}
	
	
	public OverLapCnsWithTE(String wigFile1, String wigFile2, String teFile, HashMap<String, String> maizeGenome) {
		try {
			HashMap<String, HashSet<Integer>> cnsRecords1 = new HashMap<String, HashSet<Integer>> ();
			HashMap<String, HashSet<Integer>> cnsRecords2 = new HashMap<String, HashSet<Integer>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				cnsRecords1.put(chr, new HashSet<Integer>());
				cnsRecords2.put(chr, new HashSet<Integer>());
			}
			
			readWigFile(wigFile1, cnsRecords1);
			readWigFile(wigFile2, cnsRecords2);
			HashSet<String> teClasses = new HashSet<String>();
			teClasses.add("DHH");
			teClasses.add("DTA");
			teClasses.add("DTC");
			teClasses.add("DTH");
			teClasses.add("DTM");
			teClasses.add("DTT");
			teClasses.add("DTX");
			teClasses.add("RIL");
			teClasses.add("RIT");
			teClasses.add("RLC");
			teClasses.add("RLG");
			teClasses.add("RLX");
			teClasses.add("RST");
			for( String teClass : teClasses ) {
				double totalLength = 0;
				double cns1Length = 0;
				double cns2Length = 0;
				double cns1TeLength = 0;
				double cns2TeLength = 0;
				double teLength = 0;
				
				for( String chr : cnsRecords1.keySet() ) { // so here only chromosomes 1-10 were used
	    			HashSet<Integer> teRecords = new HashSet<Integer> ();
					readGffFile(teFile, teRecords, chr, teClass);
	    			totalLength += maizeGenome.get(chr).length();
	    			teLength += teRecords.size();
	    			for( int i : cnsRecords1.get(chr) ) {
	    				if( cnsRecords2.get(chr).contains(i)  ) {
							++cns2Length;
							if( teRecords.contains(i)  ) {
		    					++cns2TeLength;
		    				}
	    				}else { // if it is not unknown, then it is known
							++cns1Length;
							if( teRecords.contains(i)  ) {
		    					++cns1TeLength;
		    				}
	    				}
	    			}
	    		}
				
				double knownCnsEWnrichmentFold = (cns1TeLength/cns1Length)/(teLength/totalLength);
	    		double unKnownCnsEWnrichmentFold = (cns2TeLength/cns2Length)/(teLength/totalLength);
	    		System.out.print(teClass +"\t" + knownCnsEWnrichmentFold +"\t" + unKnownCnsEWnrichmentFold + "\t" + cns1TeLength + "\t" + cns1Length + "\t" + teLength + "\t" + totalLength);
	    		if( knownCnsEWnrichmentFold> unKnownCnsEWnrichmentFold) {
	    			System.out.print("\tAAA");
	    		}else {
	    			System.out.print("\tBBB");
	    		}
	    		System.out.println();
				
	    		/*
				HashMap<String, Double> cns1Lengths = new HashMap<String, Double>();
				HashMap<String, Double> cns2Lengths = new HashMap<String, Double>();
				HashMap<String, Double> cns1TeLengths = new HashMap<String, Double>();
				HashMap<String, Double> cns2TeLengths = new HashMap<String, Double>();
				HashMap<String, Double> teLengths = new HashMap<String, Double>();
				
	    		for( String chr : cnsRecords1.keySet() ) { // so here only chromosomes 1-10 were used
	    			HashMap<String, HashSet<Integer>> teRecords = new HashMap<String, HashSet<Integer>> ();
					readGffFile2(teFile, teRecords, chr, teClass);
	    			for ( String family : teRecords.keySet() ) {
	    				if( !cns1Lengths.containsKey(family) ) {
	    					cns1Lengths.put(family, 0.0);
	    					cns2Lengths.put(family, 0.0);
	    					cns1TeLengths.put(family, 0.0);
	    					cns2TeLengths.put(family, 0.0);
	    					teLengths.put(family, 0.0);
	    				}
	    				teLengths.put(family, teLengths.get(family)+teRecords.get(family).size());
	    				
		    			for( int i : cnsRecords1.get(chr)) {
		    				if( cnsRecords2.get(chr).contains(i)  ) {
								cns2Lengths.put(family, cns2Lengths.get(family)+1);
								if( teRecords.get(family).contains(i)  ) {
			    					cns2TeLengths.put(family, cns2TeLengths.get(family)+1);
			    				}
		    				}else { // if it is not unknown, then it is known
								cns1Lengths.put(family, cns1Lengths.get(family)+1);
								if( teRecords.get(family).contains(i)  ) {
			    					cns1TeLengths.put(family, cns1TeLengths.get(family)+1);
			    				}
		    				}
		    			}
	    			}
	    		}
	    		
	    		
	    		for(  String family : cns1Lengths.keySet() ) {
	    			knownCnsEWnrichmentFold = (cns1TeLengths.get(family)/cns1Lengths.get(family)) / (teLengths.get(family)/totalLength);
		    		unKnownCnsEWnrichmentFold = (cns2TeLengths.get(family)/cns2Lengths.get(family)) / (teLengths.get(family)/totalLength);
		    		System.out.print("\t" + family +"\t" + knownCnsEWnrichmentFold +"\t" + unKnownCnsEWnrichmentFold + "\t" + cns1TeLengths.get(family) + "\t" + cns1Lengths.get(family) 
		    			+ "\t" + teLengths.get(family) + "\t" + totalLength);
		    		if( knownCnsEWnrichmentFold> unKnownCnsEWnrichmentFold) {
		    			System.out.print("\tAAA");
		    		}else {
		    			System.out.print("\tBBB");
		    		}
		    		System.out.println();
	    		}
	    		*/
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public OverLapCnsWithTE(String wigFile, String teFile, HashMap<String, String> maizeGenome) {
		try {
			HashMap<String, HashSet<Integer>> cnsRecords = new HashMap<String, HashSet<Integer>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				cnsRecords.put(chr, new HashSet<Integer>());
			}
			
			readWigFile(wigFile, cnsRecords);
			HashSet<String> teClasses = new HashSet<String>();
			teClasses.add("DHH");
			teClasses.add("DTA");
			teClasses.add("DTC");
			teClasses.add("DTH");
			teClasses.add("DTM");
			teClasses.add("DTT");
			teClasses.add("DTX");
			teClasses.add("RIL");
			teClasses.add("RIT");
			teClasses.add("RLC");
			teClasses.add("RLG");
			teClasses.add("RLX");
			teClasses.add("RST");
			System.out.println("teClass\tenrichmentFold\tcnsTeLength\tcnsLength\tteLength\tgenomeSize");
			for( String teClass : teClasses ) {
				double totalLength = 0;
				double cnsLength = 0;
				double cnsTeLength = 0;
				double teLength = 0;
				
				for( String chr : cnsRecords.keySet() ) { // so here only chromosomes 1-10 were used
	    			HashSet<Integer> teRecords = new HashSet<Integer> ();
					readGffFile(teFile, teRecords, chr, teClass);
	    			totalLength += maizeGenome.get(chr).length();
	    			teLength += teRecords.size();
	    			for( int i : cnsRecords.get(chr) ) {
						++cnsLength;
						if( teRecords.contains(i)  ) {
	    					++cnsTeLength;
	    				}
	    			}
	    		}
	    		double enrichmentFold = (cnsTeLength/cnsLength)/(teLength/totalLength);
	    		System.out.println(teClass +"\t" + enrichmentFold +"\t"  + cnsTeLength + "\t" + cnsLength + "\t" + teLength + "\t" + totalLength);

	    		/*
				HashMap<String, Double> cnsLengths = new HashMap<String, Double>();
				HashMap<String, Double> cnsTeLengths = new HashMap<String, Double>();
				HashMap<String, Double> teLengths = new HashMap<String, Double>();
				
	    		for( String chr : cnsRecords.keySet() ) { // so here only chromosomes 1-10 were used
	    			HashMap<String, HashSet<Integer>> teRecords = new HashMap<String, HashSet<Integer>> ();
					readGffFile2(teFile, teRecords, chr, teClass);
	    			for ( String family : teRecords.keySet() ) {
	    				if( !cnsLengths.containsKey(family) ) {
	    					cnsLengths.put(family, 0.0);
	    					cnsTeLengths.put(family, 0.0);
	    					teLengths.put(family, 0.0);
	    				}
	    				teLengths.put(family, teLengths.get(family)+teRecords.get(family).size());
	    				
		    			for( int i : cnsRecords.get(chr)) {
							cnsLengths.put(family, cnsLengths.get(family)+1);
							if( teRecords.get(family).contains(i)  ) {
		    					cnsTeLengths.put(family, cnsTeLengths.get(family)+1);
			    			}
		    			}
	    			}
	    		}
	    		
	    		for(  String family : cnsLengths.keySet() ) {
	    			enrichmentFold = (cnsTeLengths.get(family)/cnsLengths.get(family)) / (teLengths.get(family)/totalLength);
		    		System.out.print("\t" + family +"\t" + enrichmentFold + cnsTeLengths.get(family) + "\t" + cnsLengths.get(family) 
		    			+ "\t" + teLengths.get(family) + "\t" + totalLength);
		    		System.out.println();
	    		}*/
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void readGffFile(String teFile, HashSet<Integer> teRecords, String chr, String teClass) {
		try {
			File file = new File( teFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			Pattern p = Pattern.compile("(\\S+)\\s+\\S+\\s+\\S+\\s+(\\d+)\\s+(\\d+)\\s+.*ID=(\\w{3})");
			while ((tempString = reader.readLine()) != null) {
				Matcher m = p.matcher(tempString);
				if( m.find() && m.group(1).compareTo(chr)==0 && m.group(4).compareTo(teClass)==0 ) {
					for( int i=Integer.parseInt(m.group(2)); i<=Integer.parseInt(m.group(3)); ++i ) {
						teRecords.add(i);
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public void readGffFile2(String teFile, HashMap<String, HashSet<Integer>> teRecords, String chr, String teClass) {
		try {
			File file = new File( teFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			Pattern p = Pattern.compile("(\\S+)\\s+\\S+\\s+\\S+\\s+(\\d+)\\s+(\\d+)\\s+.*ID=((\\w{3})\\d+)");
			while ((tempString = reader.readLine()) != null) {
				Matcher m = p.matcher(tempString);
				if( m.find() && m.group(1).compareTo(chr)==0 && m.group(5).compareTo(teClass)==0 ) {
					if( ! teRecords.containsKey(m.group(4)) ) {
						teRecords.put(m.group(4), new HashSet<Integer>());
					}
					for( int i=Integer.parseInt(m.group(2)); i<=Integer.parseInt(m.group(3)); ++i ) {
						teRecords.get(m.group(4)).add(i);
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void readWigFile(String wigFile, HashMap<String, HashSet<Integer>> wigRecords) {
		try {
			File file = new File( wigFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			String chr = "";
			Pattern p =Pattern.compile("chrom=(\\S+)");
			while ((tempString = reader.readLine()) != null) {
				if( tempString.charAt(0) == 't' ) {
					
				}else if( tempString.charAt(0) == 'v' ) {
					Matcher m = p.matcher(tempString);
					if( m.find() ) {
						chr = m.group(1);
					}
				}else {
					String[] arrOfStr = tempString.split("\\s+");
					if( wigRecords.containsKey(chr) ) {
						if( Integer.parseInt(arrOfStr[1]) > 0 ) {
							wigRecords.get(chr).add(Integer.parseInt(arrOfStr[0]));
						}
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
