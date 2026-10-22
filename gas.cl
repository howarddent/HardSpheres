/*
  HardSpheres kernels: N hard spheres of diameter 1 and mass 1 in the box
  [0, Lx] x [0, Ly] x [0, Lz], the top face z = Lz a piston moving at W.
  One substep of length dt is

    move     x += v dt (+ any collision correction), then reflect off the walls (the piston in its own
             frame), accumulating the momentum given to the walls and the
             kinetic energy, reduced per work group into partial[group]
    clear    zero the cell counts
    bin      put each sphere's index into its cell (at most maxp)
    pairs    each sphere picks the overlapping, approaching neighbour it
             touched earliest
    apply    mutual picks collide elastically, exactly as at contact
             (impulse along the line of centres then, plus the
             displacement since, applied by the next move); the rest wait
             a substep. One collision per sphere per substep is what keeps
             energy exact: summing several impulses at once does not

  pos.w carries |v|, for colouring. pos may be a GL vertex buffer shared
  through cl_khr_gl_sharing.
*/


__kernel void move(__global float4* pos, __global float4* vel, __global float4* corr,
                   const int n, const float dt, const float4 L, const float W, const float r,
                   __global float4* partial, __local float4* scratch)
{
    int i = get_global_id(0);
    int l = get_local_id(0);
    float4 acc = (float4)(0.0f);
    if (i < n) {
        float4 p = pos[i];
        float4 v = vel[i];
        float imp = 0.0f, impP = 0.0f;
        p.xyz += v.xyz * dt + corr[i].xyz;
        corr[i] = (float4)(0.0f);
        // x and y: fixed walls
        if (p.x < r)       { p.x = 2.0f * r - p.x;       if (v.x < 0.0f) { imp -= 2.0f * v.x; v.x = -v.x; } }
        if (p.x > L.x - r) { p.x = 2.0f * (L.x - r) - p.x; if (v.x > 0.0f) { imp += 2.0f * v.x; v.x = -v.x; } }
        if (p.y < r)       { p.y = 2.0f * r - p.y;       if (v.y < 0.0f) { imp -= 2.0f * v.y; v.y = -v.y; } }
        if (p.y > L.y - r) { p.y = 2.0f * (L.y - r) - p.y; if (v.y > 0.0f) { imp += 2.0f * v.y; v.y = -v.y; } }
        // z: fixed floor at 0, the piston (the top face) at L.z moving at W
        if (p.z < r)       { p.z = 2.0f * r - p.z;       if (v.z < 0.0f) { imp -= 2.0f * v.z; v.z = -v.z; } }
        float hi = L.z - r;
        if (p.z > hi) {
            p.z = 2.0f * hi - p.z;
            if (v.z > W) {
                float nv = 2.0f * W - v.z;
                float d = v.z - nv;
                imp += d; impP += d;
                v.z = nv;
            }
        }
        p.w = length(v.xyz);
        pos[i] = p;
        vel[i] = v;
        acc = (float4)(imp, impP, 0.5f * dot(v.xyz, v.xyz), 0.0f);
    }
    scratch[l] = acc;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = get_local_size(0) / 2; s > 0; s >>= 1) {
        if (l < s) scratch[l] += scratch[l + s];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (l == 0) partial[get_group_id(0)] += scratch[0];
}

__kernel void clear(__global int* counts, const int ncells)
{
    int c = get_global_id(0);
    if (c < ncells) counts[c] = 0;
}

// cells are cubes of side cs >= 1 (the diameter), holding at most maxp
// spheres; invcs = 1/cs
inline int4 cellOf(float4 p, int4 g, float invcs)
{
    int4 c = convert_int4(floor(p * invcs));
    return clamp(c, (int4)(0), g - (int4)(1));
}

__kernel void bin(__global const float4* pos, const int n, const int4 g,
                  const float invcs, const int maxp,
                  __global int* counts, __global int* cells, __global int* overflow)
{
    int i = get_global_id(0);
    if (i >= n) return;
    int4 c = cellOf(pos[i], g, invcs);
    int cell = c.x + g.x * (c.y + g.y * c.z);
    int k = atomic_inc(&counts[cell]);
    if (k < maxp) cells[cell * maxp + k] = i;
    else atomic_inc(overflow);
}

// time since an approaching, overlapping pair touched: the root tau >= 0
// of |d - u tau| = 1, d = xj - xi, u = vj - vi (d.u < 0, |d| < 1)
inline float contactAge(float3 d, float3 u)
{
    float uu = dot(u, u);
    float du = dot(d, u);
    float disc = du * du - uu * (dot(d, d) - 1.0f);
    return (du + sqrt(max(disc, 0.0f))) / uu;
}

// each sphere's partner this substep: of the overlapping neighbours still
// approaching, the one it touched earliest (-1: none)
__kernel void pairs(__global const float4* pos, __global const float4* vel,
                    const int n, const int4 g, const float invcs, const int maxp,
                    __global const int* counts, __global const int* cells,
                    __global int* partner)
{
    int i = get_global_id(0);
    if (i >= n) return;
    float3 xi = pos[i].xyz;
    float3 vi = vel[i].xyz;
    int best = -1;
    float bestAge = -1.0f;
    int4 c = cellOf(pos[i], g, invcs);
    for (int dz = -1; dz <= 1; dz++) {
        int cz = c.z + dz;
        if (cz < 0 || cz >= g.z) continue;
        for (int dy = -1; dy <= 1; dy++) {
            int cy = c.y + dy;
            if (cy < 0 || cy >= g.y) continue;
            for (int dx = -1; dx <= 1; dx++) {
                int cx = c.x + dx;
                if (cx < 0 || cx >= g.x) continue;
                int cell = cx + g.x * (cy + g.y * cz);
                int m = min(counts[cell], maxp);
                for (int k = 0; k < m; k++) {
                    int j = cells[cell * maxp + k];
                    if (j == i) continue;
                    float3 d = pos[j].xyz - xi;
                    float r2 = dot(d, d);
                    if (r2 < 1.0f) {
                        float3 u = vel[j].xyz - vi;
                        if (dot(d, u) < 0.0f) {
                            float age = contactAge(d, u);
                            if (age > bestAge || (age == bestAge && j < best)) {
                                bestAge = age;
                                best = j;
                            }
                        }
                    }
                }
            }
        }
    }
    partner[i] = best;
}

// mutual partners collide: the elastic impulse along the line of centres
// at contact, and the displacement the new velocity would have made since
// contact (applied by the next move). Others keep their velocity and try
// again next substep. Both spheres of a pair compute mirror images, so
// momentum and energy are conserved exactly, pair by pair.
__kernel void apply(__global const float4* pos, __global const float4* vin,
                    __global float4* vout, __global float4* corr, const int n,
                    __global const int* partner, __global int* ncoll)
{
    int i = get_global_id(0);
    if (i >= n) return;
    float4 vi = vin[i];
    int j = partner[i];
    if (j >= 0 && partner[j] == i) {
        float3 d = pos[j].xyz - pos[i].xyz;
        float3 u = vin[j].xyz - vi.xyz;
        float tau = contactAge(d, u);
        float3 nc = d - u * tau;                  // unit length at contact
        nc *= rsqrt(dot(nc, nc));
        float3 dv = dot(u, nc) * nc;
        vi.xyz += dv;
        corr[i] = (float4)(dv * tau, 0.0f);
        if (i < j) atomic_inc(ncoll);
    }
    vout[i] = vi;
}

__kernel void scale(__global float4* vel, const int n, const float s)
{
    int i = get_global_id(0);
    if (i < n) vel[i] *= s;
}
