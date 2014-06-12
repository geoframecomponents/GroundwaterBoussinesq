package org.wordpress.growworkinghard.usefulClasses;

import java.text.DecimalFormat;

public class ComputePattern {

	public DecimalFormat computePattern(int endValue) {

		int[] c = new int[String.valueOf(endValue).length()];

		StringBuilder builder = new StringBuilder(c.length);

		for (int i : c) {

			builder.append(c[i]);

		}

		String pattern = builder.toString();

		return new DecimalFormat(pattern);

	}

}
