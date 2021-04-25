class PRTMaterial extends Material {
    constructor(vertexShader, fragmentShader) {

        let colorMat3 = getMat3ValueFromRGB(precomputeL[guiParams.envmapId]);

        let precomputeLR = colorMat3[0];
        let precomputeLG = colorMat3[1];
        let precomputeLB = colorMat3[2];

        super({
            'uPrecomputeLR': { type: 'precomputeL', value: precomputeLR },
            'uPrecomputeLG': { type: 'precomputeL', value: precomputeLG },
            'uPrecomputeLB': { type: 'precomputeL', value: precomputeLB },
        
        }, ['aPrecomputeLT'], vertexShader, fragmentShader, null);
    }

}

async function buildPRTMaterial(vertexPath, fragmentPath) {
    

    let vertexShader = await getShaderString(vertexPath);
    let fragmentShader = await getShaderString(fragmentPath);

    return new PRTMaterial(vertexShader, fragmentShader);

}