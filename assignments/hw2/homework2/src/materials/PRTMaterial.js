class PRTMaterial extends Material {
    constructor(vertexShader, fragmentShader) {
        let precomputeLR = mat3.create();
        let precomputeLG = mat3.create();
        let precomputeLB = mat3.create();
        

        for (let i = 0; i < 9; i++) {
            precomputeLR[i] = precomputeL[guiParams.envmapId][i][0];
            precomputeLG[i] = precomputeL[guiParams.envmapId][i][1];
            precomputeLB[i] = precomputeL[guiParams.envmapId][i][2];
        }

        super({
            'uPrecomputeLR': { type: 'matrix3fv', value: precomputeLR },
            'uPrecomputeLG': { type: 'matrix3fv', value: precomputeLG },
            'uPrecomputeLB': { type: 'matrix3fv', value: precomputeLB },
        
        }, ['aPrecomputeLT'], vertexShader, fragmentShader, null);
    }

}

async function buildPRTMaterial(vertexPath, fragmentPath) {
    

    let vertexShader = await getShaderString(vertexPath);
    let fragmentShader = await getShaderString(fragmentPath);

    return new PRTMaterial(vertexShader, fragmentShader);

}