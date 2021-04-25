function getRotationPrecomputeL(precompute_L, rotationMatrix){
	let result = [];

	let M = math.identity(9);

	// 3 by 3
	let M3 = computeSquareMatrix_3by3(rotationMatrix);
	M.subset(math.index(math.range(1, 4), math.range(1, 4)), M3);

	// 5 by 5
	let M5 = computeSquareMatrix_5by5(rotationMatrix);
	M.subset(math.index(math.range(4, 9), math.range(4, 9)), M5);

	for(var i = 0; i<3; i++) {
		let L = [precompute_L[0][i], precompute_L[1][i], precompute_L[2][i],
				 precompute_L[3][i], precompute_L[4][i], precompute_L[5][i],
				 precompute_L[6][i], precompute_L[7][i], precompute_L[8][i]];
		
		let vecL = math.matrix(L);

		let ML = math.multiply(M, vecL);
		ML = ML.reshape([9]);
		
		result[i] = [ML.subset(math.index(0)), ML.subset(math.index(1)), ML.subset(math.index(2)),
					ML.subset(math.index(3)), ML.subset(math.index(4)), ML.subset(math.index(5)),
					ML.subset(math.index(6)), ML.subset(math.index(7)), ML.subset(math.index(8)),];
	}

	return result;
}

function computeSquareMatrix_3by3(rotationMatrix){ // 计算方阵SA(-1) 3*3 
	
	// 1、pick ni - {ni}
	let n1 = [1, 0, 0, 0]; let n2 = [0, 0, 1, 0]; let n3 = [0, 1, 0, 0];

	// 2、{P(ni)} - A  A_inverse
	let N = math.zeros(3, 3);
	for (let i = 0; i < 3; i++) {
		N.subset(math.index(i, 0), n1[i]);
		N.subset(math.index(i, 1), n2[i]);
		N.subset(math.index(i, 2), n3[i]);
	}

	// project vector to SH
	let P1 = SHEval(n1[0], n1[1], n1[2], 3);
	let P2 = SHEval(n2[0], n2[1], n2[2], 3);
	let P3 = SHEval(n3[0], n3[1], n3[2], 3);

	let A = math.zeros(3, 3);
	for (let i = 0; i < 3; i++) {
		A.subset(math.index(i, 0), P1[i+1]);
		A.subset(math.index(i, 1), P2[i+1]);
		A.subset(math.index(i, 2), P3[i+1]);
	}

	let A_inv = math.inv(A);

	// 3、用 R 旋转 ni - {R(ni)}
	let R = mat4Matrix2mathMatrix(rotationMatrix);
	R = R.subset(math.index(math.range(0, 3), math.range(0, 3)));

	let Rn= math.multiply(math.transpose(R), N);

	// 4、R(ni) SH投影 - S
	// project rotated vector to SH
	let PR1 = SHEval(Rn.subset(math.index(0, 0)), Rn.subset(math.index(1, 0)), Rn.subset(math.index(2, 0)), 3);
	let PR2 = SHEval(Rn.subset(math.index(0, 1)), Rn.subset(math.index(1, 1)), Rn.subset(math.index(2, 1)), 3);
	let PR3 = SHEval(Rn.subset(math.index(0, 2)), Rn.subset(math.index(1, 2)), Rn.subset(math.index(2, 2)), 3);

	let S = math.zeros(3, 3);
	for (let i = 0; i < 3; i++) {
		S.subset(math.index(i, 0), PR1[i+1]);
		S.subset(math.index(i, 1), PR2[i+1]);
		S.subset(math.index(i, 2), PR3[i+1]);
	}

	// 5、S*A_inverse
	let M = math.multiply(S, A_inv);

	return M;
}

function computeSquareMatrix_5by5(rotationMatrix){ // 计算方阵SA(-1) 5*5
	
	// 1、pick ni - {ni}
	let k = 1 / math.sqrt(2);
	let n1 = [1, 0, 0, 0]; let n2 = [0, 0, 1, 0]; let n3 = [k, k, 0, 0]; 
	let n4 = [k, 0, k, 0]; let n5 = [0, k, k, 0];

	// 2、{P(ni)} - A  A_inverse
	let N = math.zeros(3, 5);
	for (let i = 0; i < 3; i++) {
		N.subset(math.index(i, 0), n1[i]);
		N.subset(math.index(i, 1), n2[i]);
		N.subset(math.index(i, 2), n3[i]);
		N.subset(math.index(i, 3), n4[i]);
		N.subset(math.index(i, 4), n5[i]);
	}

	// project vector to SH
	let P1 = SHEval(n1[0], n1[1], n1[2], 3);
	let P2 = SHEval(n2[0], n2[1], n2[2], 3);
	let P3 = SHEval(n3[0], n3[1], n3[2], 3);
	let P4 = SHEval(n4[0], n4[1], n4[2], 3);
	let P5 = SHEval(n5[0], n5[1], n5[2], 3);

	let A = math.zeros(5, 5);
	for (let i = 0; i < 5; i++) {
		A.subset(math.index(i, 0), P1[i+4]);
		A.subset(math.index(i, 1), P2[i+4]);
		A.subset(math.index(i, 2), P3[i+4]);
		A.subset(math.index(i, 3), P4[i+4]);
		A.subset(math.index(i, 4), P5[i+4]);
	}

	let A_inv = math.inv(A);

	// 3、用 R 旋转 ni - {R(ni)}
	let R = mat4Matrix2mathMatrix(rotationMatrix);
	R = R.subset(math.index(math.range(0, 3), math.range(0, 3)));

	let Rn= math.multiply(math.transpose(R), N);

	// 4、R(ni) SH投影 - S
	let PR1 = SHEval(Rn.subset(math.index(0, 0)), Rn.subset(math.index(1, 0)), Rn.subset(math.index(2, 0)), 3);
	let PR2 = SHEval(Rn.subset(math.index(0, 1)), Rn.subset(math.index(1, 1)), Rn.subset(math.index(2, 1)), 3);
	let PR3 = SHEval(Rn.subset(math.index(0, 2)), Rn.subset(math.index(1, 2)), Rn.subset(math.index(2, 2)), 3);
	let PR4 = SHEval(Rn.subset(math.index(0, 3)), Rn.subset(math.index(1, 3)), Rn.subset(math.index(2, 3)), 3);
	let PR5 = SHEval(Rn.subset(math.index(0, 4)), Rn.subset(math.index(1, 4)), Rn.subset(math.index(2, 4)), 3);

	let S = math.zeros(5, 5);
	for (let i = 0; i < 5; i++) {
		S.subset(math.index(i, 0), PR1[i+4]);
		S.subset(math.index(i, 1), PR2[i+4]);
		S.subset(math.index(i, 2), PR3[i+4]);
		S.subset(math.index(i, 3), PR4[i+4]);
		S.subset(math.index(i, 4), PR5[i+4]);
	}

	// 5、S*A_inverse
	let M = math.multiply(S, A_inv);

	return M;
}

function mat4Matrix2mathMatrix(rotationMatrix){

	let mathMatrix = [];
	for(let i = 0; i < 4; i++){
		let r = [];
		for(let j = 0; j < 4; j++){
			r.push(rotationMatrix[i*4+j]);
		}
		mathMatrix.push(r);
	}
	return math.matrix(mathMatrix);

}

function getMat3ValueFromRGB(precomputeL){

    let colorMat3 = [];
    for(var i = 0; i<3; i++){
        colorMat3[i] = mat3.fromValues( precomputeL[0][i], precomputeL[1][i], precomputeL[2][i],
										precomputeL[3][i], precomputeL[4][i], precomputeL[5][i],
										precomputeL[6][i], precomputeL[7][i], precomputeL[8][i] ); 
	}
    return colorMat3;
}