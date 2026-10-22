/*
  StickySpheres kernels: N spheres of diameter 1 and mass 1 in the box
  [0, Lx] x [0, Ly] x [0, Lz] (the top, z = Lz, a piston moving at W
  along z; gravity g along -z), interacting through a stepped potential:

    a hard core at r = R[0] = 1, then K flat shells, shell s spanning
    R[s] <= r < R[s+1] at energy eps U[s]; zero beyond the cutoff R[K]

  A square well is K = 1; a stepped Lennard-Jones potential is K ~ 9,
  each U[s] the LJ energy averaged over the shell (see ustickysim.pas).
  Every force is an impulse at a step: within a shell spheres move freely.

  Each sphere keeps a list of its neighbours within the cutoff and the
  shell each is in (at most MAXB, set at build time; the first nn[i]
  slots are used, kept compact). The lists of a pair always agree.

  Pair events, along the line of centres n at the moment they happen,
  with u_n the relative normal velocity and u_n^2/4 the pair's relative
  kinetic energy along n (reduced mass 1/2), and dU the potential energy
  the step would add:

    core     r reaches 1 (shell 0), approaching: elastic bounce
    inward   r reaches R[s] from shell s (s > 0), or the cutoff R[K]
             from outside (the pair joins the lists), approaching:
               capture on, one of the pair free (no bonds), the step
               downhill (dU < 0) and u_n^2/4 < eps:
                 the pair steps in with u_n unchanged and -dU is given off
                 (inelastic: a slow pair falls into the well and sticks);
               otherwise, if u_n^2/4 > dU, it steps in, u_n^2/4 changing
                 by -dU; if not (an uphill step it can't climb), it
                 bounces back
    outward  r reaches R[s+1] from shell s, separating: if u_n^2/4 > dU it
             steps out (beyond the cutoff: leaves the lists), otherwise
             bounces back

  A "bond" is a pair in a shell below BondEdge (r < R[BondEdge], 1.5 for
  Lennard-Jones): those are counted, coloured, drawn and make clusters.

  One substep of length dt:
    move   x += v dt (+ the correction of the last event), gravity, walls;
           colour (pos.w = bonds + min(|v|/vscale, 0.999)); reduce per
           work group: wall impulse, piston impulse, KE, bond ends, sum of
           heights, pair potential energy (half per list entry)
    clear, bin    the cell list (cells at least R[K] wide)
    pairs  each sphere picks its earliest pending event
    apply  mutual picks happen, exactly as at the event; each sphere
           updates only its own velocity and list, so nothing races
    pairs and apply are repeated (sweeps) within a substep until no event
    is pending (after the first sweep, only for the spheres that had one or
    were picked as a partner: the rest can't have gained one) (pairs counts the spheres with one): a sphere in a cluster
    often has several due at once - in the stepped Lennard-Jones well a
    bound pair rattles in shells 0.06 wide - and a fixed number of sweeps
    let a backlog snowball in dense clusters, pairs drifting far outside
    the shells they were listed in.
  The floor can attract: a one-body well of depth epsF for centres below
  zw (one diameter above contact), so a liquid can wet it and spread.
  Gravity accelerates both spheres of a pair alike, so their relative
  motion - all the pair events use - is a straight line.
*/

#ifndef MAXB
#define MAXB 48
#endif

__kernel void move(__global float4* pos, __global float4* vel, __global float4* corr,
                   const int n, const float dt, const float4 L, const float W, const float r,
                   __global const int* nn, __global const int* shl,
                   __global const float* U, const int bondEdge, const float eps,
                   const float vscale, const float g, const float epsF, const float zw,
                   __global float8* partial, __local float8* scratch)
{
    int i = get_global_id(0);
    int l = get_local_id(0);
    float8 acc = (float8)(0.0f);
    if (i < n) {
        float4 p = pos[i];
        float4 v = vel[i];
        float imp = 0.0f, impP = 0.0f;
        float oldz = p.z;
        p.xyz += v.xyz * dt + corr[i].xyz;
        p.z -= 0.5f * g * dt * dt;
        v.z -= g * dt;
        corr[i] = (float4)(0.0f);
        if (p.x < r)       { p.x = 2.0f * r - p.x;       if (v.x < 0.0f) { imp -= 2.0f * v.x; v.x = -v.x; } }
        if (p.x > L.x - r) { p.x = 2.0f * (L.x - r) - p.x; if (v.x > 0.0f) { imp += 2.0f * v.x; v.x = -v.x; } }
        if (p.y < r)       { p.y = 2.0f * r - p.y;       if (v.y < 0.0f) { imp -= 2.0f * v.y; v.y = -v.y; } }
        if (p.y > L.y - r) { p.y = 2.0f * (L.y - r) - p.y; if (v.y > 0.0f) { imp += 2.0f * v.y; v.y = -v.y; } }
        // z: the floor at 0, the piston (the top face) at L.z moving at W
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
        // the floor's attraction: a well of depth epsF for centres below zw.
        // Crossing zw downwards a sphere gains epsF of energy along z;
        // upwards it loses it, or bounces back if it hasn't enough. tau is
        // how long ago it crossed (the position after it follows the new
        // velocity), capped as the pair corrections are. If the move
        // crossed against v.z (a pair correction did it) the sphere is put
        // back on its side instead
        // The well's pull is a force between floor and sphere too, so its
        // momentum changes count in the wall impulse (negative: towards
        // the floor); without them a populated well made the floor's
        // pressure read many times too high
        float inWell = 0.0f;
        if (epsF > 0.0f) {
            float vz0 = v.z;
            int was = oldz < zw, now = p.z < zw;
            if (was && !now) {
                if (v.z <= 0.0f) p.z = 2.0f * zw - p.z;
                else if (v.z * v.z > 2.0f * epsF) {
                    float tau = min((p.z - zw) / v.z, 2.0f * dt);
                    float nv = sqrt(v.z * v.z - 2.0f * epsF);
                    p.z = zw + nv * tau;
                    v.z = nv;
                } else {
                    p.z = 2.0f * zw - p.z;
                    v.z = -v.z;
                }
            } else if (!was && now) {
                if (v.z >= 0.0f) p.z = 2.0f * zw - p.z;
                else {
                    float tau = min((p.z - zw) / v.z, 2.0f * dt);
                    float nv = -sqrt(v.z * v.z + 2.0f * epsF);
                    p.z = max(zw + nv * tau, r);
                    v.z = nv;
                }
            }
            imp += v.z - vz0;
            if (p.z < zw) inWell = 1.0f;
        }
        int m = nn[i];
        int nbond = 0;
        float pe = 0.0f;
        for (int k = 0; k < m; k++) {
            int s = shl[i * MAXB + k];
            pe += U[s];
            if (s < bondEdge) nbond++;
        }
        p.w = (float)nbond + min(length(v.xyz) / vscale, 0.999f);
        pos[i] = p;
        vel[i] = v;
        acc = (float8)(imp, impP, 0.5f * dot(v.xyz, v.xyz), (float)nbond, p.z, 0.5f * eps * pe, inWell, 0.0f);
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

// time since an approaching pair, now inside radius sqrt(R2), crossed it
inline float rootIn(float3 d, float3 u, float R2)
{
    float uu = dot(u, u), du = dot(d, u);
    float disc = du * du - uu * (dot(d, d) - R2);
    return (du + sqrt(max(disc, 0.0f))) / uu;
}

// where and when a separating pair, now at |d| > R, crossed R: tau ago,
// nc the unit line of centres then. From the straight-line history when
// that crossed R within the last 0.05 diameters of relative travel;
// otherwise (a history bent by other events, or no crossing on it) along
// the current line of centres, tau = (|d| - R)/u_n. The closest-approach
// point once used there gave u_n = 0 and pairs drifted apart still listed;
// the current line every time overestimates u_n and lets pairs escape too
// easily
inline float exitGeom(float3 d, float3 u, float R, float3* nc)
{
    float uu = dot(u, u), du = dot(d, u);
    float disc = du * du - uu * (dot(d, d) - R * R);
    if (disc > 0.0f) {
        float tau = (du - sqrt(disc)) / uu;
        if (tau * tau * uu < 0.0025f) {
            *nc = (d - u * tau) / R;
            return tau;
        }
    }
    float r = length(d);
    *nc = d / r;
    return (r - R) / max(du / r, 1e-6f);
}

// event types
#define EV_CORE 1
#define EV_ENTER 2     // reaches the cutoff from outside
#define EV_IN 3        // a step inward, between shells
#define EV_OUT 4       // a step outward (from the last shell: leaves)

// the pending event of a pair and how long ago it happened. s: its shell
// (-1: not listed). Symmetric in the pair (d, u -> -d, -u), so both
// spheres find the same event
inline int eventOf(float3 d, float3 u, int s, int canAdd, int K,
                   __global const float* R, float* age)
{
    float dd = dot(d, d), du = dot(d, u);
    if (s >= 0) {
        float r0 = R[s];
        if (dd < r0 * r0 && du < 0.0f) {
            *age = rootIn(d, u, r0 * r0);
            return s == 0 ? EV_CORE : EV_IN;
        }
        float r1 = R[s + 1];
        if (dd > r1 * r1 && du > 0.0f) {
            float3 nc;
            *age = exitGeom(d, u, r1, &nc);
            return EV_OUT;
        }
        return 0;
    }
    float rc = R[K];
    if (canAdd && dd < rc * rc && du < 0.0f) { *age = rootIn(d, u, rc * rc); return EV_ENTER; }
    if (dd < 1.0f && du < 0.0f) { *age = rootIn(d, u, 1.0f); return EV_CORE; }
    return 0;
}

// slot of j in i's list, or -1
inline int slotOf(__global const int* nbr, __global const int* nn, int i, int j)
{
    int m = nn[i];
    for (int k = 0; k < m; k++)
        if (nbr[i * MAXB + k] == j) return k;
    return -1;
}

__kernel void pairs(__global const float4* pos, __global const float4* vel,
                    const int n, const int4 g, const float invcs, const int maxp,
                    __global const int* counts, __global const int* cells,
                    __global const int* nbr, __global const int* shl, __global const int* nn,
                    __global const float* R, const int K,
                    __global int* partner, __global int* ptype, __global int* pending,
                    __global const int* active, __global int* activeNext)
{
    int i = get_global_id(0);
    if (i >= n) return;
    // after the first sweep only spheres that had an event pending, or were
    // picked as a partner, can have one: the rest skip the search
    if (!active[i]) {
        partner[i] = -1;
        ptype[i] = 0;
        return;
    }
    float3 xi = pos[i].xyz;
    float3 vi = vel[i].xyz;
    int best = -1, bestType = 0;
    float bestAge = -1.0f, age;
    // listed neighbours: steps between shells, the core, leaving
    int m = nn[i];
    for (int k = 0; k < m; k++) {
        int j = nbr[i * MAXB + k];
        int t = eventOf(pos[j].xyz - xi, vel[j].xyz - vi, shl[i * MAXB + k], 0, K, R, &age);
        if (t && (age > bestAge || (age == bestAge && j < best))) {
            bestAge = age; best = j; bestType = t;
        }
    }
    // the rest: reaching the cutoff (or, lists full, the core)
    float rc2 = R[K] * R[K];
    int canI = m < MAXB;
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
                int mc = min(counts[cell], maxp);
                for (int k = 0; k < mc; k++) {
                    int j = cells[cell * maxp + k];
                    if (j == i) continue;
                    float3 d = pos[j].xyz - xi;
                    if (dot(d, d) >= rc2) continue;
                    float3 u = vel[j].xyz - vi;
                    if (dot(d, u) >= 0.0f) continue;      // only approaching pairs can enter
                    if (slotOf(nbr, nn, i, j) >= 0) continue;
                    int t = eventOf(d, u, -1, canI && nn[j] < MAXB, K, R, &age);
                    if (t && (age > bestAge || (age == bestAge && j < best))) {
                        bestAge = age; best = j; bestType = t;
                    }
                }
            }
        }
    }
    partner[i] = best;
    ptype[i] = bestType;
    if (best >= 0) {
        atomic_inc(pending);
        activeNext[i] = 1;
        activeNext[best] = 1;
    }
}

__kernel void fill(__global int* a, const int n, const int v)
{
    int i = get_global_id(0);
    if (i < n) a[i] = v;
}

// mutual picks happen. cnt: [0] core bounces, [1] bonds formed (steps in
// across R[BondEdge]), [2] bonds broken, [3] captures, [4 + b] captures at
// boundary R[b] (for the energy given off); each pair counted once
__kernel void apply(__global const float4* pos, __global const float4* vin,
                    __global float4* vout, __global float4* corr, const int n,
                    __global const int* partner, __global const int* ptype,
                    __global int* nbr, __global int* shl, __global int* nn,
                    __global const float* R, __global const float* U, const int K,
                    const int bondEdge, const float eps, const int capture,
                    const float taumax, __global int* cnt)
{
    int i = get_global_id(0);
    if (i >= n) return;
    float4 vi = vin[i];
    int j = partner[i];
    int type = ptype[i];
    if (j >= 0 && partner[j] == i && ptype[j] == type) {
        float3 d = pos[j].xyz - pos[i].xyz;
        float3 u = vin[j].xyz - vi.xyz;
        int first = i < j;
        int slot = (type == EV_ENTER) ? -1 : slotOf(nbr, nn, i, j);
        int s = (slot >= 0) ? shl[i * MAXB + slot] : K;   // K: outside
        // the boundary crossed, the shell beyond it, and the energy step
        int b, snew;
        float tau;
        float3 nc;
        if (type == EV_CORE) {
            b = 0;
            snew = s;
            tau = rootIn(d, u, 1.0f);
            nc = d - u * tau;
            nc *= rsqrt(dot(nc, nc));
        } else if (type == EV_OUT) {
            b = s + 1;
            snew = s + 1;
            tau = exitGeom(d, u, R[b], &nc);
        } else {
            b = (type == EV_ENTER) ? K : s;
            snew = b - 1;
            tau = rootIn(d, u, R[b] * R[b]);
            nc = d - u * tau;
            nc *= rsqrt(dot(nc, nc));
        }
        float un = dot(u, nc);
        float un2 = un * un;
        float unew;
        int moved = 0;
        if (type == EV_CORE) {
            unew = -un;
            if (first) atomic_inc(&cnt[0]);
        } else {
            float uold = (s < K) ? U[s] : 0.0f;
            float unw = (snew < K) ? U[snew] : 0.0f;
            float dU = eps * (unw - uold);
            // capture stands in for the third body that carries a gas
            // pair's energy away, so it applies only when one of the two
            // is free (no bonds as of this substep's move: pos.w < 1). In
            // a cluster the neighbours do that by ordinary collisions;
            // capturing there too drained the liquid's heat (the
            // thermostat's uniform rescaling put it into the vapour) and
            // froze the liquid solid
            int gas = pos[i].w < 1.0f || pos[j].w < 1.0f;
            if (type != EV_OUT && capture && gas && dU < 0.0f && 0.25f * un2 < eps) {
                unew = un;                    // steps in, giving off -dU
                moved = 1;
                if (first) { atomic_inc(&cnt[3]); atomic_inc(&cnt[4 + b]); }
            } else if (0.25f * un2 > dU) {
                float w = sqrt(max(un2 - 4.0f * dU, 0.0f));
                unew = (type == EV_OUT) ? w : -w;
                moved = 1;
            } else
                unew = -un;                   // can't climb the step: bounces
        }
        if (moved) {
            if (snew >= K) {
                // leaves the cutoff: off the list, the last entry into its slot
                int last = nn[i] - 1;
                nbr[i * MAXB + slot] = nbr[i * MAXB + last];
                shl[i * MAXB + slot] = shl[i * MAXB + last];
                nbr[i * MAXB + last] = -1;
                nn[i] = last;
            } else if (slot < 0) {
                int k = nn[i];
                nbr[i * MAXB + k] = j;
                shl[i * MAXB + k] = snew;
                nn[i] = k + 1;
            } else
                shl[i * MAXB + slot] = snew;
            if (first && b == bondEdge) atomic_inc(&cnt[(type == EV_OUT) ? 2 : 1]);
        }
        float3 dv = -0.5f * (unew - un) * nc;
        vi.xyz += dv;
        // the displacement the new velocity would have made since the
        // event, tau capped at two substeps (a barely separating pair gives
        // a huge tau, and the correction threw spheres across the box)
        corr[i] += (float4)(dv * min(tau, taumax), 0.0f);
    }
    vout[i] = vi;
}

__kernel void scale(__global float4* vel, const int n, const float s)
{
    int i = get_global_id(0);
    if (i < n) vel[i] *= s;
}
