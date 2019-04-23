R"(

float GetEVDWIJ(float rij, float epsij, float roij) {
	float r6ij = pow((roij / rij), 6);
	return epsij * (pow(r6ij, 2) - 2.0 * r6ij);
}

float GetEElstIJ(float rij, float qi, float qj, float epsilon, float ceu_to_kcal) {
	return ceu_to_kcal * qi * qj / (epsilon * rij);
}

kernel void CalculateNonBondedForces(constant int *combinationMatrix, constant int *nonInts, constant int *nonIntsSize, constant float3 *positions, constant float *charges, constant float *vdwRadii, constant float *vdwAttractionMagnitudes, constant float *dielectric, constant float *CEU_TO_KCAL, global float *eVDW, global float *eElst) {
	int id = get_global_id(0);

	int i = combinationMatrix[id * 2];
	int j = combinationMatrix[id * 2 + 1];

	int found = 0;

	for (int k = 0; k < *nonIntsSize; k += 2) {
		if (nonInts[k] == i && nonInts[k + 1] == j) {
			found = 1;
		}
	}

	if (!found) {
		float a = positions[i].x - positions[j].x;
		float b = positions[i].y - positions[j].y;
		float c = positions[i].z - positions[j].z;

		float distance = sqrt(a * a + b * b + c * c);
		float vdwAttractionMagnitudeIJ = sqrt(vdwAttractionMagnitudes[i]) * sqrt(vdwAttractionMagnitudes[j]);
		float vdwRadiusIJ = vdwRadii[i] + vdwRadii[j];

		(*eVDW) += GetEVDWIJ(distance, vdwAttractionMagnitudeIJ, vdwRadiusIJ);
		(*eElst) += GetEElstIJ(distance, charges[i], charges[j], *dielectric, *CEU_TO_KCAL);
	}
}
)"