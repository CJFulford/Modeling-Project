
// output to be interpolated between vertices and passed to the fragment stage
out vec2 textureCoords;

//uniforms
// add transformation uniforms
uniform vector sphereVec;
uniform vector triangleVec;
uniform vector planeVec;

out vec3 vertexPoint;
out vec3 center;

void main()
{
    // assign vertex position without modification
    gl_Position = proj * view * model * vec4(VertexPosition, 1.0);
	center = vec3(model[3][0], model[3][1], model[3][2]);
	
    vec4 p = model * vec4(VertexPosition, 1.0);
	vertexPoint = vec3(p.x, p.y, p.z);

    textureCoords = VertexTextureCoords;
}
