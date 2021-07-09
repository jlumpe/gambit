"""GAMBIT 0.1.0

Revision ID: c43540b80d50
Revises: 
Create Date: 2021-07-08 13:34:30.131392

Creates 0.1.0 database from scratch.
"""
from alembic import op
import sqlalchemy as sa

from gambit.db.sqla import JsonString


# revision identifiers, used by Alembic.
revision = 'c43540b80d50'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.create_table('genome_sets',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('key', sa.String(), nullable=False),
        sa.Column('version', sa.String(), nullable=True),
        sa.Column('name', sa.String(), nullable=False),
        sa.Column('description', sa.String(), nullable=True),
        sa.Column('extra', JsonString(), nullable=True),
        sa.PrimaryKeyConstraint('id', name=op.f('pk_genome_sets')),
        sa.UniqueConstraint('key', 'version', name=op.f('uq_genome_sets_key'))
    )
    op.create_index(op.f('ix_genome_sets_key'), 'genome_sets', ['key'], unique=False)

    op.create_table('genomes',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('key', sa.String(), nullable=False),
        sa.Column('description', sa.String(), nullable=False),
        sa.Column('ncbi_db', sa.String(), nullable=True),
        sa.Column('ncbi_id', sa.Integer(), nullable=True),
        sa.Column('genbank_acc', sa.String(), nullable=True),
        sa.Column('refseq_acc', sa.String(), nullable=True),
        sa.Column('extra', JsonString(), nullable=True),
        sa.PrimaryKeyConstraint('id', name=op.f('pk_genomes')),
        sa.UniqueConstraint('genbank_acc', name=op.f('uq_genomes_genbank_acc')),
        sa.UniqueConstraint('key', name=op.f('uq_genomes_key')),
        sa.UniqueConstraint('ncbi_db', 'ncbi_id', name=op.f('uq_genomes_ncbi_db')),
        sa.UniqueConstraint('refseq_acc', name=op.f('uq_genomes_refseq_acc'))
    )

    op.create_table('taxa',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('key', sa.String(), nullable=False),
        sa.Column('name', sa.String(), nullable=False),
        sa.Column('rank', sa.String(), nullable=True),
        sa.Column('description', sa.String(), nullable=True),
        sa.Column('distance_threshold', sa.Float(), nullable=True),
        sa.Column('report', sa.Boolean(), server_default=sa.text('1'), nullable=False),
        sa.Column('genome_set_id', sa.Integer(), nullable=False),
        sa.Column('parent_id', sa.Integer(), nullable=True),
        sa.Column('ncbi_id', sa.Integer(), nullable=True),
        sa.Column('extra', JsonString(), nullable=True),
        sa.ForeignKeyConstraint(['genome_set_id'], ['genome_sets.id'], name=op.f('fk_taxa_genome_set_id_genome_sets'), ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['parent_id'], ['taxa.id'], name=op.f('fk_taxa_parent_id_taxa'), ondelete='SET NULL'),
        sa.PrimaryKeyConstraint('id', name=op.f('pk_taxa')),
        sa.UniqueConstraint('key', name=op.f('uq_taxa_key'))
    )
    op.create_index(op.f('ix_taxa_genome_set_id'), 'taxa', ['genome_set_id'], unique=False)
    op.create_index(op.f('ix_taxa_name'), 'taxa', ['name'], unique=False)
    op.create_index(op.f('ix_taxa_ncbi_id'), 'taxa', ['ncbi_id'], unique=False)
    op.create_index(op.f('ix_taxa_parent_id'), 'taxa', ['parent_id'], unique=False)
    op.create_index(op.f('ix_taxa_rank'), 'taxa', ['rank'], unique=False)

    op.create_table('genome_annotations',
        sa.Column('genome_id', sa.Integer(), nullable=False),
        sa.Column('genome_set_id', sa.Integer(), nullable=False),
        sa.Column('taxon_id', sa.Integer(), nullable=True),
        sa.Column('organism', sa.String(), nullable=True),
        sa.ForeignKeyConstraint(['genome_id'], ['genomes.id'], name=op.f('fk_genome_annotations_genome_id_genomes'), ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['genome_set_id'], ['genome_sets.id'], name=op.f('fk_genome_annotations_genome_set_id_genome_sets'), ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['taxon_id'], ['taxa.id'], name=op.f('fk_genome_annotations_taxon_id_taxa'), ondelete='SET NULL'),
        sa.PrimaryKeyConstraint('genome_id', 'genome_set_id', name=op.f('pk_genome_annotations'))
    )
    op.create_index(op.f('ix_genome_annotations_taxon_id'), 'genome_annotations', ['taxon_id'], unique=False)


def downgrade():
    op.drop_index(op.f('ix_genome_annotations_taxon_id'), table_name='genome_annotations')
    op.drop_table('genome_annotations')
    op.drop_index(op.f('ix_taxa_rank'), table_name='taxa')
    op.drop_index(op.f('ix_taxa_parent_id'), table_name='taxa')
    op.drop_index(op.f('ix_taxa_ncbi_id'), table_name='taxa')
    op.drop_index(op.f('ix_taxa_name'), table_name='taxa')
    op.drop_index(op.f('ix_taxa_genome_set_id'), table_name='taxa')
    op.drop_table('taxa')
    op.drop_table('genomes')
    op.drop_index(op.f('ix_genome_sets_key'), table_name='genome_sets')
    op.drop_table('genome_sets')
