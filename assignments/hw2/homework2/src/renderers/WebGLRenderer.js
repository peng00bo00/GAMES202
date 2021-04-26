class WebGLRenderer {
    meshes = [];
    shadowMeshes = [];
    lights = [];

    constructor(gl, camera) {
        this.gl = gl;
        this.camera = camera;
    }

    addLight(light) {
        this.lights.push({
            entity: light,
            meshRender: new MeshRender(this.gl, light.mesh, light.mat)
        });
    }
    addMeshRender(mesh) { this.meshes.push(mesh); }
    addShadowMeshRender(mesh) { this.shadowMeshes.push(mesh); }

    render() {
        const gl = this.gl;

        gl.clearColor(0.0, 0.0, 0.0, 1.0); // Clear to black, fully opaque
        gl.clearDepth(1.0); // Clear everything
        gl.enable(gl.DEPTH_TEST); // Enable depth testing
        gl.depthFunc(gl.LEQUAL); // Near things obscure far things

        console.assert(this.lights.length != 0, "No light");
        console.assert(this.lights.length == 1, "Multiple lights");

        const timer = Date.now() * 0.0005;

        for (let l = 0; l < this.lights.length; l++) {
            // Draw light
            this.lights[l].meshRender.mesh.transform.translate = this.lights[l].entity.lightPos;
            this.lights[l].meshRender.draw(this.camera);

            // Shadow pass
            if (this.lights[l].entity.hasShadowMap == true) {
                for (let i = 0; i < this.shadowMeshes.length; i++) {
                    this.shadowMeshes[i].draw(this.camera);
                }
            }

            // Camera pass
            for (let i = 0; i < this.meshes.length; i++) {
                this.gl.useProgram(this.meshes[i].shader.program.glShaderProgram);
                this.gl.uniform3fv(this.meshes[i].shader.program.uniforms.uLightPos, this.lights[l].entity.lightPos);

                for (let k in this.meshes[i].material.uniforms) {

                    let cameraModelMatrix = mat4.create();
                    // mat4.fromRotation(cameraModelMatrix, timer, [0, 1, 0]);

                    if (k == 'uMoveWithCamera') { // The rotation of the skybox
                        gl.uniformMatrix4fv(
                            this.meshes[i].shader.program.uniforms[k],
                            false,
                            cameraModelMatrix);
                    }

                    // Bonus - Fast Spherical Harmonic Rotation
                    let precomputeL_RGBMat3 = getRotationPrecomputeL(precomputeL[guiParams.envmapId], cameraModelMatrix);

                    if (k == 'uPrecomputeLR') {
                        let precomputeLR = mat3.fromValues(precomputeL_RGBMat3[0][0], precomputeL_RGBMat3[0][1], precomputeL_RGBMat3[0][2],
                                                            precomputeL_RGBMat3[0][3], precomputeL_RGBMat3[0][4], precomputeL_RGBMat3[0][5],
                                                            precomputeL_RGBMat3[0][6], precomputeL_RGBMat3[0][7], precomputeL_RGBMat3[0][8]);

                        gl.uniformMatrix3fv(
                            this.meshes[i].shader.program.uniforms[k],
                            false,
                            precomputeLR);
                    }

                    if (k == 'uPrecomputeLG') {
                        let precomputeLG = mat3.fromValues(precomputeL_RGBMat3[1][0], precomputeL_RGBMat3[1][1], precomputeL_RGBMat3[1][2],
                                                            precomputeL_RGBMat3[1][3], precomputeL_RGBMat3[1][4], precomputeL_RGBMat3[1][5],
                                                            precomputeL_RGBMat3[1][6], precomputeL_RGBMat3[1][7], precomputeL_RGBMat3[1][8]);

                        gl.uniformMatrix3fv(
                            this.meshes[i].shader.program.uniforms[k],
                            false,
                            precomputeLG);
                    }

                    if (k == 'uPrecomputeLB') {
                        let precomputeLB = mat3.fromValues(precomputeL_RGBMat3[2][0], precomputeL_RGBMat3[2][1], precomputeL_RGBMat3[2][2],
                                                            precomputeL_RGBMat3[2][3], precomputeL_RGBMat3[2][4], precomputeL_RGBMat3[2][5],
                                                            precomputeL_RGBMat3[2][6], precomputeL_RGBMat3[2][7], precomputeL_RGBMat3[2][8]);

                        gl.uniformMatrix3fv(
                            this.meshes[i].shader.program.uniforms[k],
                            false,
                            precomputeLB);
                    }
                    
                }

                this.meshes[i].draw(this.camera);
            }
        }

    }
}