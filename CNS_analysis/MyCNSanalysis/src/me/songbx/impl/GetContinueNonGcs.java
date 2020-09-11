package me.songbx.impl;

import java.util.ArrayList;

import me.songbx.model.ContinueNonGc;

public class GetContinueNonGcs {
	
	public static ArrayList<ContinueNonGc> getContinueNonGcs( String seq ){
		ArrayList<ContinueNonGc> continueNonGcs = new ArrayList<ContinueNonGc>();
		int position=1;
		boolean lastNonGc = false;
		int length = 0;
		seq = seq.toUpperCase();
		for( int i=0; i<seq.length(); ++i ) {
			if( seq.charAt(i) == 'G' || seq.charAt(i) == 'C' ) {
				if( lastNonGc ) {
					ContinueNonGc continueNonGc = new ContinueNonGc(position, length);
					continueNonGcs.add(continueNonGc);
				}
				lastNonGc = false;
				length = 0;
			} else {
				if( lastNonGc ) {
					++length;
				}else {
					length = 1;
					position = i + 1;
				}
				lastNonGc = true;
			}
		}
		if( length > 0 ) {
			ContinueNonGc continueNonGc = new ContinueNonGc(position, length);
			continueNonGcs.add(continueNonGc);
		}
		
		return continueNonGcs;
	}

	public static void main(String[] args) {
		ArrayList<ContinueNonGc> continueNonGcs = GetContinueNonGcs.getContinueNonGcs("ATCGATGCATGCATGC");
		for ( ContinueNonGc continueNonGc : continueNonGcs ){
			System.out.println(continueNonGc.getPosition() + "\t" + continueNonGc.getLength());
		}
	}
}
