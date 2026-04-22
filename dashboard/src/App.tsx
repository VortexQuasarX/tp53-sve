import { useState } from 'react'
import './index.css'
import { Login } from './components/Login'
import { Overview } from './components/Overview'
import { Leaderboard } from './components/Leaderboard'

type TabId = 'dashboard' | 'resume' | 'interviews' | 'leaderboard' | 'retention'

const TABS: { id: TabId; label: string; icon: string }[] = [
  { id: 'dashboard', label: 'Triage Overview', icon: '🧬' },
  { id: 'structure', label: '3D Structure View', icon: '🔬' },
  { id: 'kan', label: 'KAN Explainability', icon: '🧠' },
  { id: 'leaderboard', label: 'Variant Ranking', icon: '🏆' },
  { id: 'vaccine', label: 'mRNA Pipeline', icon: '🧪' },
]

export default function App() {
  const [user, setUser] = useState<any>(null)
  const [tab, setTab] = useState<TabId>('dashboard')

  if (!user) {
    return <Login onLogin={setUser} />
  }

  return (
    <div className="portal-container">
      {/* Sidebar */}
      <div className="sidebar">
        <div className="sidebar-logo">
          <div className="logo-icon">p53</div>
          <span>TP53-SVE</span>
        </div>
        
        <div className="nav-list">
          {TABS.map(t => (
            <button 
              key={t.id} 
              className={`nav-item ${tab === t.id ? 'active' : ''}`} 
              onClick={() => setTab(t.id)}
            >
              <span className="nav-icon">{t.icon}</span>
              <span className="nav-label">{t.label}</span>
            </button>
          ))}
        </div>

        <div className="sidebar-footer">
          <div className="user-pill">
            <img src={user.avatar} alt="avatar" />
            <div className="user-info">
              <span className="user-name">{user.name}</span>
              <span className="user-role">Student</span>
            </div>
          </div>
        </div>
      </div>

      {/* Main Content */}
      <div className="main-viewport">
        <header className="main-header">
          <div className="breadcrumb">TP53-SVE / {TABS.find(t => t.id === tab)?.label}</div>
          <button className="primary-btn">Rerun Analysis</button>
        </header>

        <div className="content-area">
          {tab === 'dashboard' && <Overview />}
          {tab === 'leaderboard' && <Leaderboard />}
          {/* Other tabs placeholder */}
          {['structure', 'kan', 'vaccine'].includes(tab) && (
            <div className="placeholder-view">
              <h2>{TABS.find(t => t.id === tab)?.label}</h2>
              <p>Scientific module integration in progress...</p>
            </div>
          )}
        </div>
      </div>
    </div>
  )
}
