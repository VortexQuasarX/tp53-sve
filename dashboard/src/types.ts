export interface Mutation {
    m: string;
    rmsd: number;
    cls: string;
    crit: string;
    pos: number;
    wt: string;
    mut: string;
}

export interface SummaryData {
    mean: number;
    max: number;
    maxRes: number;
    above5: number;
    above10: number;
}

export interface ToolComparison {
    m: string;
    rmsd: number;
    sev: string;
    sift: number;
    pp2: number;
    rRank: number;
    sRank: number;
    pRank: number;
}

export interface DomainData {
    m: string;
    full: number;
    nterm: number;
    dbd: number;
    tetra: number;
    cterm: number;
}

export interface CorrelationData {
    m: string;
    rmsd: number;
    plddt: number;
    charge: number;
    mass: number;
    hydro: number;
    dbdDist: number;
}

export type PerResidueMap = Record<string, number[]>;
export type SummaryMap = Record<string, SummaryData>;

export type TabId = 'overview' | 'ranking' | 'residue' | 'criteria' | 'compare' | 'position' | 'tools' | 'insights';

export const SEVERITY_COLORS: Record<string, string> = {
    Critical: '#ef4444',
    Severe: '#f97316',
    Moderate: '#eab308',
    Stable: '#22c55e',
};

export const CRITERION_COLORS: Record<string, string> = {
    Phase1: '#3b82f6',
    A: '#a855f7',
    B: '#22c55e',
    C: '#f97316',
    D: '#ef4444',
    E: '#eab308',
    F: '#94a3b8',
};

export const CRITERION_LABELS: Record<string, string> = {
    Phase1: 'Phase 1 Hotspot',
    A: 'Same-Position',
    B: 'Benign Control',
    C: 'Non-DBD Domain',
    D: 'Gain-of-Function',
    E: 'Temp-Sensitive',
    F: 'Rare Pathogenic',
};

export function getSeverity(rmsd: number): string {
    if (rmsd > 30) return 'Critical';
    if (rmsd > 20) return 'Severe';
    if (rmsd > 10) return 'Moderate';
    return 'Stable';
}
