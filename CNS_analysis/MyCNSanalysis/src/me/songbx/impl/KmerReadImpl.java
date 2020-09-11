package me.songbx.impl;

import java.io.*;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author song
 * @version 1.0, 2019-12-07
 */
public class KmerReadImpl {
	// k-mer_seq, frequency in the genome/sequence, GC content
	private HashMap<String, Integer> kmers = new HashMap<String, Integer>();
	public KmerReadImpl(String fastaFileLocation) {
		
		File file = new File(fastaFileLocation);
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
						kmers.put(seq, Integer.parseInt(sequenceName) );
					}
					sequenceName = m1.group(1);
					sequenceBuffer = new StringBuffer();
				} else {
					sequenceBuffer.append(tempString);
				}
			}
			if (sequenceName.length() > 0 && sequenceBuffer.length() > 0) {
				String seq = sequenceBuffer.toString().toUpperCase();
				seq=seq.replaceAll("\\s", "");
				kmers.put(seq, Integer.parseInt(sequenceName) );
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
		}
	}
	public HashMap<String, Integer> getKmers() {
		return kmers;
	}
	public void setKmers(HashMap<String, Integer> kmers) {
		this.kmers = kmers;
	}
	public static double getGcContent(String seq) {
		String s = seq.toUpperCase();
		double gc=0.0;
		double seqLength = 0.0;
		for( int i=0; i<s.length(); ++i ) {
			if( s.charAt(i) == 'C' || s.charAt(i) == 'G' ) {
				gc += 1.0;
				seqLength += 1.0;
			} else if( s.charAt(i) == 'A' || s.charAt(i) == 'T' ) {
				seqLength += 1.0;
			}
		}
		return gc/seqLength; 
	}
}
