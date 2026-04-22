import React from 'react'

export const Leaderboard: React.FC = () => {
  const mutations = [
    { rank: 1, name: 'R175H', score: 99.8, trend: 'up' },
    { rank: 2, name: 'R248Q', score: 98.2, trend: 'stable' },
    { rank: 3, name: 'G245S', score: 95.4, trend: 'up' },
    { rank: 4, name: 'R273H', score: 92.1, trend: 'down' },
    { rank: 5, name: 'Y220C', score: 88.5, trend: 'up' },
  ]

  return (
    <div className="leaderboard-module">
      <div className="module-header">
        <h1>Pathogenicity Ranking</h1>
        <p>A comprehensive ranking of clinical TP53 variants based on KAN-predicted structural impact.</p>
      </div>

      <div className="card">
        <div className="tbl-wrap">
          <table>
            <thead>
              <tr>
                <th>Rank</th>
                <th>Mutation Variant</th>
                <th>KAN Impact %</th>
                <th>Trend</th>
              </tr>
            </thead>
            <tbody>
              {mutations.map(m => (
                <tr key={m.rank} className={m.name === 'R175H' ? 'my-rank' : ''}>
                  <td>
                    <div className={`rank-circle rank-${m.rank}`}>
                      {m.rank}
                    </div>
                  </td>
                  <td>
                    <div className="player-cell">
                      <div className="avatar-small">🧬</div>
                      <span>{m.name} {m.name === 'R175H' && '(Target)'}</span>
                    </div>
                  </td>
                  <td>
                    <div className="score-cell">
                      <b>{m.score}</b>
                    </div>
                  </td>
                  <td>
                    <span className={`trend-icon ${m.trend}`}>
                      {m.trend === 'up' ? '▲' : m.trend === 'down' ? '▼' : '▬'}
                    </span>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      <style>{`
        .rank-circle {
          width: 32px; height: 32px;
          border-radius: 50%;
          display: grid; place-items: center;
          font-weight: 800; font-size: 0.8rem;
          background: #f3f4f6; color: var(--text-grey);
        }
        .rank-1 { background: #fef3c7; color: #92400e; border: 1px solid #fcd34d; }
        .rank-2 { background: #e0e7ff; color: #3730a3; border: 1px solid #c7d2fe; }
        .rank-3 { background: #ffedd5; color: #9a3412; border: 1px solid #fed7aa; }

        .player-cell { display: flex; align-items: center; gap: 12px; }
        .avatar-small { 
          width: 28px; height: 28px; 
          border-radius: 8px; 
          background: var(--accent); 
          color: white; 
          display: grid; place-items: center; 
          font-size: 0.75rem; font-weight: 700;
        }
        
        .my-rank { background: var(--accent-light) !important; }
        .my-rank td { border-bottom-color: var(--accent) !important; }

        .trend-icon.up { color: #10b981; }
        .trend-icon.down { color: #ef4444; }
        .trend-icon.stable { color: #9ca3af; }

        th { background: #f9fafb; color: var(--text-grey); font-size: 0.7rem; text-transform: uppercase; padding: 12px 16px; }
        td { padding: 12px 16px; border-bottom: 1px solid var(--border-light); font-size: 0.85rem; }
      `}</style>
    </div>
  )
}
