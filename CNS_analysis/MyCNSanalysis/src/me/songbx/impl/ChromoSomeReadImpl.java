package me.songbx.impl;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import me.songbx.model.FastaIndexEntry;
import me.songbx.model.Strand;

/**
 * @author song
 * @version 1.0, 2014-07-09
 */
public class ChromoSomeReadImpl {
	/**hashmap of chromosomes, the index are the names of chromosomes, the value are the chromosomes*/
	private HashMap<String, String> chromoSomeHashMap = new HashMap<String, String>();
	public ChromoSomeReadImpl() {
		
	}
	/**
	 * @param the location of fasta or multiple fasta format file. The meta-line should be like ">Chr2"
	 */
	public ChromoSomeReadImpl(String fastaFileLocation) {
		FastaIndexEntryImpl fastaIndexEntryImpl = new FastaIndexEntryImpl();
		RandomAccessFile reader = null;

		String fastaIndexFileLocation = fastaFileLocation + ".fai";
		File f = new File(fastaIndexFileLocation);
		if(f.exists() && !f.isDirectory()) {
			fastaIndexEntryImpl.readFastaIndexFile(fastaIndexFileLocation);
		}else {
			fastaIndexEntryImpl.createFastaIndexFile(fastaFileLocation);
		}
		File file = new File(fastaFileLocation);
		try {
			reader = new RandomAccessFile(file, "r");

			for ( String name : fastaIndexEntryImpl.getEntries().keySet() ) {
				FastaIndexEntry entry = fastaIndexEntryImpl.getEntries().get(name);
				int newlines_in_sequence = entry.getLength()/entry.getLine_blen();
				int seqlen = newlines_in_sequence + entry.getLength();
				byte[] b = new byte[seqlen];
				try {
					reader.seek(entry.getOffset());
					reader.read(b);
				} catch (IOException e) {
					e.printStackTrace();
				}
				String seq = new String(b);
				seq = seq.replaceAll("\\s", "");

				chromoSomeHashMap.put(name, seq);
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


		/*
		File file = new File(fileLocation);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			Pattern p1 = Pattern.compile("^>(\\S+)");
			StringBuffer sequenceBuffer = new StringBuffer();
			String sequenceName = "";
			while ((tempString = reader.readLine()) != null) {
				
				Matcher m1 = p1.matcher(tempString);
				if (m1.find()) {
					if (sequenceName.length() > 0 && sequenceBuffer.length() > 0) {
						String seq = sequenceBuffer.toString().toUpperCase();
						seq=seq.replaceAll("\\s", "");
						ChromoSome c = new ChromoSome(sequenceName, seq);
						chromoSomeHashMap.put(sequenceName, c);
					}
					sequenceName = m1.group(1);
					sequenceBuffer = new StringBuffer();
					//System.out.println(sequenceName);
				} else {
					sequenceBuffer.append(tempString);
				}
			}
			if (sequenceName.length() > 0 && sequenceBuffer.length() > 0) {
				String seq = sequenceBuffer.toString().toUpperCase();
				seq=seq.replaceAll("\\s", "");
				ChromoSome c = new ChromoSome(sequenceName, seq);
				chromoSomeHashMap.put(sequenceName, c);
				//System.out.println(sequenceName);
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e1) {

				}
			}
		}*/
	}
	
	/**
	 * get a chromosome with a name, if could find this name, return the chromosome, if could not return null
	 * @param ChromoSome name
	 * @return ChromoSome
	 */
	public synchronized String getChromoSomeById(String name) {
		if (chromoSomeHashMap.containsKey(name)) {
			return chromoSomeHashMap.get(name);
		}
		return null;
	}

	/**
	 * if start1 > end1, switch them.
	 * And then, the sequence position counts start from 1, the start1-th character would be returned, the end1-th would not be returned
	 * please make sure start >= 1, end <= the length of chromosome sequence
	 * 
	 * @param chromeSomeName
	 * @param start1
	 *            : the start position of positive sequence 
	 * @param end1
	 *            : the start position of positive sequence
	 * @param strand
	 *            : want positive(POSITIVE) sequence, or reverse complementary
	 *            sequence (NEGATIVE)
	 * @return the wanted sequence
	 */
	public synchronized String getSubSequence(String chromeSomeName, int start1, int end1,
			Strand strand) {
		int start;
		int end;
		if (start1 < end1) {
			start = start1;
			end = end1;
		} else {
			start = end1;
			end = start1;
		}

		String subSequence = "";
		if (chromoSomeHashMap.containsKey(chromeSomeName)) {
			String cS = chromoSomeHashMap.get(chromeSomeName);
			start--;
			if(start<0){
				start=0;
			}
			if( end > cS.length()){
				end = cS.length();
			}
			if( start > cS.length()){
				start = cS.length();
			}
			if( end<0 ){
				end=0;
			}
			//System.out.println("start: " + start + " end: " + end);
			subSequence = cS.subSequence(start, end).toString();
			if (Strand.NEGTIVE == strand) {
				return getReversecomplementary(subSequence);
			}
		}
		return subSequence;
	}
	
	
	
	public synchronized char getChar(String chromeSomeName, int position) {
		position--;
		String cS = chromoSomeHashMap.get(chromeSomeName);
		char subSequence = cS.charAt(position);
		return subSequence;
	}

	/**
	 * @param DNA
	 *            /RNA sequence
	 * @return reverse complementary DNA sequence(only T, no U)
	 */
	public synchronized static String getReversecomplementary(String sequence) {
		StringBuffer reversecomplementary = new StringBuffer();
		for (int i = sequence.length() - 1; i >= 0; i--) {
			char c = sequence.charAt(i);
			if ('A' == c) {
				c = 'T';
			} else if ('T' == c) {
				c = 'A';
			} else if ('U' == c) {
				c = 'A';
			} else if ('C' == c) {
				c = 'G';
			} else if ('G' == c) {
				c = 'C';
			} else if ('R' == c) {
				c = 'Y';
			} else if ('Y' == c) {
				c = 'R';
			} else if ('K' == c) {
				c = 'M';
			} else if ('M' == c) {
				c = 'K';
			} else if ('B' == c) {
				c = 'V';
			} else if ('V' == c) {
				c = 'B';
			} else if ('D' == c) {
				c = 'H';
			} else if ('H' == c) {
				c = 'D';
			}
			reversecomplementary.append(c);
		}
		return reversecomplementary.toString();
	}

	public synchronized HashMap<String, String> getChromoSomeHashMap() {
		return chromoSomeHashMap;
	}

	public synchronized void setChromoSomeHashMap(
			HashMap<String, String> chromoSomeHashMap) {
		this.chromoSomeHashMap = chromoSomeHashMap;
	}
	public void releaseRam() {
		ArrayList<String> keys = new ArrayList<String>();
		Iterator<String> i = chromoSomeHashMap.keySet().iterator();
		while( i.hasNext() ){
			String name = i.next();
			chromoSomeHashMap.put(name, null);
			keys.add(name);
		}
		for(String n : keys){
			chromoSomeHashMap.remove(n);
		}
		chromoSomeHashMap=null;
	}
	
	public static void main( String[] argv ){
		ChromoSomeReadImpl chromoSomeReadImpl = new ChromoSomeReadImpl("F:\\chi_polish\\chi_v1k.fa");
		System.out.println(chromoSomeReadImpl.getChromoSomeById("NSCAFB.163"));
		String t = chromoSomeReadImpl.getSubSequence("NSCAFB.163", 200,  269, Strand.POSITIVE);
		System.out.println(t);
	}
}
