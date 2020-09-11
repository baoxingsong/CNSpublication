package me.songbx.model;

public class BedRecord implements Comparable<BedRecord>{
	int start;
	int end;
	Strand strand;
	public Strand getStrand() {
		return strand;
	}
	public void setStrand(Strand strand) {
		this.strand = strand;
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
	public BedRecord(int start, int end) {
		super();
		this.start = start;
		this.end = end;
	}
	public BedRecord(int start, int end, Strand strand) {
		super();
		this.start = start;
		this.end = end;
		this.strand=strand;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + start;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		BedRecord other = (BedRecord) obj;
		if (end != other.end)
			return false;
		if (start != other.start)
			return false;
		return true;
	}
	@Override
	public int compareTo(BedRecord arg0) {
		return this.getStart() - arg0.getStart();
	}
}
