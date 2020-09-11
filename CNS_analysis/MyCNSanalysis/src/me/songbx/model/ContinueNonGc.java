package me.songbx.model;

public class ContinueNonGc {
	private int position;
	private int length;
	public int getPosition() {
		return position;
	}
	public void setPosition(int position) {
		this.position = position;
	}
	public int getLength() {
		return length;
	}
	public void setLength(int length) {
		this.length = length;
	}
	public ContinueNonGc(int position, int length) {
		super();
		this.position = position;
		this.length = length;
	}
}
