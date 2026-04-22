import React from 'react'

export const Overview: React.FC = () => {
  const kpis = [
    { label: 'Batch Rank', value: '1629' },
    { label: 'SASA Change', value: '+420' },
    { label: 'Kan Probability', value: '0.94' },
    { label: 'Variants Analyzed', value: '128' },
    { label: 'Mean Pathogenicity', value: '0.62' },
    { label: 'Model Confidence', value: '0.98' },
  ]

  return (
    <div className="dashboard-overview">
      <div className="welcome-section">
        <div className="welcome-text">
          <h1>Welcome, Visva R</h1>
          <p>Saturday, 11 April 2026</p>
        </div>
      </div>

      <div className="kpi-grid">
        {kpis.map((kpi, i) => (
          <div key={i} className="kpi-card">
            <span className="kpi-label">{kpi.label}</span>
            <span className="kpi-value">{kpi.value}</span>
          </div>
        ))}
      </div>

      <div className="g2">
        <div className="card">
          <h3><span className="ico">🧬</span> Biological Insights</h3>
          <ul className="improvement-list">
            <li>Verify Zinc-binding pocket displacement stability.</li>
            <li>Refine SASA exposure thresholds for non-DBD variants.</li>
            <li>Analyze allosteric rewiring in the tetramer interface.</li>
          </ul>
        </div>
        
        <div className="card">
          <h3><span className="ico">📈</span> Technical Score</h3>
          <div className="chart-placeholder">
            <div className="score-ring">
              <span className="score-num">45</span>
              <span className="score-total">/100</span>
            </div>
          </div>
        </div>
      </div>

      <style>{`
        .welcome-section { margin-bottom: 32px; }
        .welcome-text h1 { font-size: 2rem; font-weight: 700; color: var(--text-main); margin-bottom: 4px; }
        .welcome-text p { color: var(--text-grey); font-size: 0.95rem; }
        
        .improvement-list { list-style: none; padding: 0; }
        .improvement-list li { 
          padding: 12px 0; 
          border-bottom: 1px solid var(--border-light); 
          font-size: 0.85rem; 
          color: var(--text-grey);
          display: flex;
          align-items: center;
          gap: 10px;
        }
        .improvement-list li::before { content: '•'; color: var(--accent); font-weight: bold; }
        .improvement-list li:last-child { border: none; }

        .chart-placeholder { height: 160px; display: grid; place-items: center; }
        .score-ring {
          width: 120px; height: 120px;
          border-radius: 50%;
          border: 8px solid #f3f4f6;
          border-top-color: var(--accent);
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
        }
        .score-num { font-size: 1.8rem; font-weight: 800; color: var(--accent); line-height: 1; }
        .score-total { font-size: 0.75rem; color: var(--text-grey); font-weight: 600; }
      `}</style>
    </div>
  )
}
