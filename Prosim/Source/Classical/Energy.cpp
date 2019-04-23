#include "Energy.h"

#include <iterator>
#include <math.h>

#include "Constants.h"
#include "Geometry.h"

#include "Utils/IterationTools.h"

namespace classical {

	static const char *s_nonBondedForcesKernelSource = {
#include "../../Kernels/NonBondedForces.cl"
	};

	NonBondedEnergyGPUCalculator::NonBondedEnergyGPUCalculator() {
		cl_uint numPlatforms; //the NO. of platforms
		cl_platform_id platform = NULL; //the chosen platform
		cl_int status = clGetPlatformIDs(0, NULL, &numPlatforms);
		if (status != CL_SUCCESS)
		{
			std::cout << "Non-bonded energy GPU calculator error: Getting platforms!" << std::endl;
			/* return FAILURE; */
		}

		/*For clarity, choose the first available platform. */
		if (numPlatforms > 0)
		{
			cl_platform_id* platforms =
				(cl_platform_id*)malloc(numPlatforms * sizeof(cl_platform_id));
			status = clGetPlatformIDs(numPlatforms, platforms, NULL);
			platform = platforms[0];
			free(platforms);
		}

		cl_uint numDevices = 0;
		cl_device_id        *devices;
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
		if (numDevices == 0) //no GPU available.
		{
			std::cout << "No GPU device available for non-bonded calculations." << std::endl;
			std::cout << "Choose CPU as default device." << std::endl;
			status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numDevices);
			devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
			status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices, devices, NULL);
		}
		else
		{
			devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
			status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
		}


		context = clCreateContext(NULL, 1, devices, NULL, NULL, NULL);

		commandQueue = clCreateCommandQueue(context, devices[0], 0, NULL);

		size_t sourceSize[] = { strlen(s_nonBondedForcesKernelSource) };
		program = clCreateProgramWithSource(context, 1, &s_nonBondedForcesKernelSource, sourceSize, NULL);

		status = clBuildProgram(program, 1, devices, NULL, NULL, NULL);

		if (status != CL_SUCCESS)
		{
			size_t len;

			printf("Error: Failed to build program executable for non-bonded calculations!\n");
			clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
			printf("Got 1st len = %lu\n", len);
			char* buffer = (char*)malloc(len);
			clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_LOG, len, buffer, &len);
			printf("Got 2nd len = %lu\n", len);
			printf("</OpenCL compiler error message>\n%s\n</OpenCL compiler error message>\n", buffer);

		}

		kernel = clCreateKernel(program, "CalculateNonBondedForces", NULL);
	}

	NonBondedEnergyGPUCalculator::~NonBondedEnergyGPUCalculator() {
		clReleaseKernel(kernel); //Release kernel.
		clReleaseProgram(program); //Release the program object.
		clReleaseMemObject(kCombinationMatrixMem); //Release mem object.
		clReleaseMemObject(kNonIntsMem);
		clReleaseMemObject(kNonIntsSizeMem);
		clReleaseMemObject(kPositionsMem);
		clReleaseMemObject(kChargesMem);
		clReleaseMemObject(kVdwRadiiMem);
		clReleaseMemObject(kVdwAttractionMagnitudesMem);
		clReleaseMemObject(kDielectricMem);
		clReleaseMemObject(kCeuToKCalMem);
		clReleaseMemObject(keElstMem);
		clReleaseMemObject(keVDWMem);
		clReleaseCommandQueue(commandQueue); //Release  Command queue.
		clReleaseContext(context); //Release context.
	}

	math::Vec2 NonBondedEnergyGPUCalculator::GetENonBonded(const std::vector<Atom *> &atoms, std::vector<int> &nonInts, double dielectric) {
		float eVDW = 0.0;
		float eElst = 0.0;

		int natoms = atoms.size();

		utils::IterationMatrix matrix(utils::CombinationsNR(natoms, 2), 2);

		utils::CombinationKN(matrix, 2, natoms);

		int kNonIntsSize = nonInts.size();

		std::vector<math::Vec3> kPositions;
		std::vector<float> kCharges;
		std::vector<float> kVdwRadii;
		std::vector<float> kVdwAttractionMagnitudes;

		for (int n = 0; n < natoms; n++) {
			kPositions.push_back(atoms[n]->position);
			kCharges.push_back((float)atoms[n]->charge);
			kVdwRadii.push_back((float)atoms[n]->vdwRadius);
			kVdwAttractionMagnitudes.push_back((float)atoms[n]->vdwAttractionMagnitude);
		}

		float kCeuToKCal = CEU_TO_KCAL;

		kCombinationMatrixMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(matrix.GetDataWritable().size()) * sizeof(int), (void *)&matrix.GetDataWritable()[0], NULL);
		kNonIntsMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(nonInts.size()) * sizeof(int), (void *)&nonInts[0], NULL);
		kNonIntsSizeMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(sizeof(int)), &kNonIntsSize, NULL);
		kPositionsMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(natoms * sizeof(math::Vec3)), &kPositions[0], NULL);
		kChargesMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(natoms * sizeof(float)), &kCharges[0], NULL);
		kVdwRadiiMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(natoms * sizeof(float)), &kVdwRadii[0], NULL);
		kVdwAttractionMagnitudesMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(natoms * sizeof(float)), &kVdwAttractionMagnitudes[0], NULL);
		kDielectricMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(sizeof(float)), &dielectric, NULL);
		kCeuToKCalMem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			(sizeof(float)), &kCeuToKCal, NULL);
		keVDWMem = clCreateBuffer(context, CL_MEM_READ_WRITE,
			sizeof(float), NULL, NULL);
		keElstMem = clCreateBuffer(context, CL_MEM_READ_WRITE,
			sizeof(float), NULL, NULL);

		clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&kCombinationMatrixMem);
		clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&kNonIntsMem);
		clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&kNonIntsSizeMem);
		clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&kPositionsMem);
		clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&kChargesMem);
		clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&kVdwRadiiMem);
		clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&kVdwAttractionMagnitudesMem);
		clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&kDielectricMem);
		clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&kCeuToKCalMem);
		clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&keVDWMem);
		clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&keElstMem);

		size_t global_work_size[1] = { matrix.GetRows() / 2 };
		clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL,
			global_work_size, NULL, 0, NULL, NULL);

		clEnqueueReadBuffer(commandQueue, keVDWMem, CL_TRUE, 0,
			sizeof(float), &eVDW, 0, NULL, NULL);
		clEnqueueReadBuffer(commandQueue, keElstMem, CL_TRUE, 0,
			sizeof(float), &eElst, 0, NULL, NULL);

		return math::Vec2(eVDW, eElst);
	}

	double GetEBond(double rij, double req, double kb) {
		return kb * pow((rij - req), 2);
	}

	double GetEAngle(double aijk, double aeq, double ka) {
		return ka * pow((DEGREES_TO_RADIANS * (aijk - aeq)), 2);
	}

	double GetETorsion(double tijkl, double vn, double gamma, int nfold, int paths) {
		return vn * (1.0 + cos(DEGREES_TO_RADIANS * (nfold * tijkl - gamma))) / paths;
	}

	double GetEOutOfPlane(double oijkl, double vn) {
		return vn * (1.0 + cos(DEGREES_TO_RADIANS * (2.0 * oijkl - 180.00)));
	}

	double GetEVDWIJ(double rij, double epsij, double roij) {
		double r6ij = pow((roij / rij), 6);
		return epsij * (pow(r6ij, 2) - 2.0 * r6ij);
	}

	double GetEElstIJ(double rij, double qi, double qj, double epsilon) {
		return CEU_TO_KCAL * qi * qj / (epsilon * rij);
	}

	double GetEBoundI(double kBox, double bound, const math::Vec3 &position, const math::Vec3 &origin, const String &boundType) {
		double eBoundI = 0.0;

		if (boundType == "cube") {
			for (int j = 0; j < 3; j++) {
				double scale = (float)(abs(position[j] - origin[j]) >= bound);
				eBoundI += scale * kBox * pow((abs(position[j] - origin[j]) - bound), 2);
			}
		} else if (boundType == "sphere") {
			double rio = GetRij(origin, position);
			double scale = (float)(rio >= bound);
			eBoundI += scale * kBox * pow((rio - bound), 2);
		}

		return eBoundI;
	}

	double GetEKineticI(double mass, const math::Vec3 &velocity) {
		double eKinI = 0.0;
		eKinI += mass * pow(velocity.x, 2);
		eKinI += mass * pow(velocity.y, 2);
		eKinI += mass * pow(velocity.z, 2);
		return 0.5 * KINETIC_TO_KCAL * eKinI;
	}

	double GetEBonds(const std::vector<Bond *> &bonds) {
		double eBonds = 0.0;

		for (Bond *bond : bonds) {
			bond->CalculateEnergy();
			eBonds += bond->energy;
		}

		return eBonds;
	}

	double GetEAngles(const std::vector<Angle *> &angles) {
		double eAngles = 0.0;

		for (Angle *angle : angles) {
			angle->CalculateEnergy();
			eAngles += angle->energy;
		}

		return eAngles;
	}

	double GetETorsions(const std::vector<Torsion *> &torsions) {
		double eTorsions = 0.0;

		for (Torsion *torsion : torsions) {
			torsion->CalculateEnergy();
			eTorsions += torsion->energy;
		}

		return eTorsions;
	}

	double GetEOutOfPlanes(const std::vector<OutOfPlane *> &outOfPlanes) {
		double eOutOfPlanes = 0.0;

		for (OutOfPlane *outOfPlane : outOfPlanes) {
			outOfPlane->CalculateEnergy();
			eOutOfPlanes += outOfPlane->energy;
		}

		return eOutOfPlanes;
	}

	math::Vec2 GetENonBonded(const std::vector<Atom *> &atoms, std::vector<int> &nonInts, double dielectric) {
		float eVDW = 0.0;
		float eElst = 0.0;

		int natoms = atoms.size();

		utils::IterationMatrix matrix(utils::CombinationsNR(natoms, 2), 2);

		utils::CombinationKN(matrix, 2, natoms);

#ifdef PS_OPTIMIZED
#pragma omp parallel for
#endif
		for (int n = 0; n < utils::CombinationsNR(natoms, 2); n++) {
			int i = matrix(n, 0);
			int j = matrix(n, 1);

			for (int k = 0; k < nonInts.size(); k += 2) {
				if (nonInts[k] == i && nonInts[k + 1] == j) {
					continue;
				}
			}

			Atom *atom1 = atoms[i];
			Atom *atom2 = atoms[j];

			double distance = GetRij(atom1->position, atom2->position);
			double vdwAttractionMagnitudeIJ = sqrt(atom1->vdwAttractionMagnitude) * sqrt(atom2->vdwAttractionMagnitude);
			double vdwRadiusIJ = atom1->vdwRadius + atom2->vdwRadius;
			eElst += GetEElstIJ(distance, atom1->charge, atom2->charge, dielectric);
			eVDW += GetEVDWIJ(distance, vdwAttractionMagnitudeIJ, vdwRadiusIJ);
		}

		return math::Vec2(eVDW, eElst);
	}

	double GetEBound(const std::vector<Atom *> &atoms, double kBox, double boundary, const math::Vec3 &origin, const String &boundType) {
		double eBound = 0.0;

		for (Atom *atom : atoms) {
			eBound += GetEBoundI(kBox, boundary, atom->position, origin, boundType);
		}

		return eBound;
	}

	double GetEKinetic(const std::vector<Atom *> &atoms, const String &kineticType) {
		if (kineticType == "noKinetic") {
			return 0.0;
		}

		double eKinetic = 0.0;
		if (kineticType == "leapfrog") {
			for (Atom *atom : atoms) {
				math::Vec3 velocity(atom->velocity.x + atom->previousVelocity.x, 
					atom->velocity.y + atom->previousVelocity.y,
					atom->velocity.z + atom->previousVelocity.z);

				velocity.Multiply(0.5);

				eKinetic += GetEKineticI(atom->mass, velocity);
			}
		}
		else {
			for (Atom *atom : atoms) {
				eKinetic += GetEKineticI(atom->mass, atom->velocity);
			}
		}

		return eKinetic;
	}

	double GetTemperature(double eKinetic, int natoms) {
		return (2.0 / 3.0) * eKinetic / (BOLTZMANN_CONSTANT * natoms);
	}

}